#include "betaEventAction.hh"

#include "G4AnalysisManager.hh"
#include "G4Event.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include "G4SDManager.hh"
#include "G4GenericMessenger.hh"
#include "G4HCofThisEvent.hh"

#include "CalorimeterHit.hh"
#include "TargetHit.hh"

// #include "betaLineage.hh"

#include <algorithm>
#include <string>
#include <sstream>
#include <iostream>

betaEventAction::betaEventAction()
    : G4UserEventAction()
{
  fDE_MeV = new std::vector<double>(225, 0.0);
  fPID = new std::vector<int>(225, 0);
  fMessenger = std::make_unique<G4GenericMessenger>(this, "/beta/", "beta controls");

  fMessenger->DeclareProperty("doTargetTrace", fDoTargetTrace,
                              "Print target dE breakdown every event (0/1).");

  fMessenger->DeclareProperty("saveRndmOnPeak", fSaveRndmOnPeak,
                              "Save RNG engine status when targetE is in [peakEmin, peakEmax] MeV (0/1).");

  fMessenger->DeclareProperty("peakEmin", fPeakEmin_MeV,
                              "Peak window min (MeV).");

  fMessenger->DeclareProperty("peakEmax", fPeakEmax_MeV,
                              "Peak window max (MeV).");

  fMessenger->DeclareProperty("rndmDir", fRndmDir,
                              "Directory to save rndm status files.");

  fMessenger->DeclareProperty("rndmPrefix", fRndmPrefix,
                              "Prefix of rndm status file name.");
}

betaEventAction::~betaEventAction() {
    // G4cout << "DEBUG: ~betaEventAction" << G4endl;
}

void betaEventAction::BeginOfEventAction(const G4Event *evt)
{
  const int eid = evt->GetEventID();

  // if (fWatchEventID >= 0 && eid == fWatchEventID) {
  //   const std::string fname = "rndm_begin_evt" + std::to_string(eid) + ".rndm";
  //   G4Random::saveEngineStatus(fname.c_str());
  //   G4cout << "[RNDM-SAVE] " << fname << "\n";
  // }

  // betaLineage::Instance()->Clear();

  fEvtEdepCell = 0.;
  fEvtEdepTarget = 0.;

  std::fill(fDE_MeV->begin(), fDE_MeV->end(), 0.0);
  std::fill(fPID->begin(), fPID->end(), 0);
}

void betaEventAction::EndOfEventAction(const G4Event *evt)
{
  auto *ana = G4AnalysisManager::Instance();
  const auto eventID = evt->GetEventID();

  auto *hce = evt->GetHCofThisEvent();

  // -------------------------
  // Calorimeter (ntupleId = 1)
  // -------------------------
  if (hce)
  {
    static G4int calHCID = -1;
    if (calHCID < 0)
    {
      calHCID = G4SDManager::GetSDMpointer()->GetCollectionID("CalorimeterHC");
    }

    auto *calHC = static_cast<CalorimeterHitsCollection *>(hce->GetHC(calHCID));
    if (calHC)
    {
      for (size_t i = 0; i < calHC->GetSize(); i++)
      {
        const auto *hit = (*calHC)[i];

        // per-hit row
        ana->FillNtupleIColumn(1, 0, eventID);
        ana->FillNtupleIColumn(1, 1, hit->GetCopyNo());
        ana->FillNtupleDColumn(1, 2, hit->GetTime() / ns);
        ana->FillNtupleDColumn(1, 3, hit->GetEdep() / MeV);
        ana->FillNtupleIColumn(1, 4, hit->GetPDG());
        ana->FillNtupleDColumn(1, 5, hit->GetMomentum().x() / MeV);
        ana->FillNtupleDColumn(1, 6, hit->GetMomentum().y() / MeV);
        ana->FillNtupleDColumn(1, 7, hit->GetMomentum().z() / MeV);
        ana->FillNtupleSColumn(1, 8, hit->GetCreator());
        ana->FillNtupleIColumn(1, 9, hit->GetOriginType());
        ana->AddNtupleRow(1);

        // calarr vector (copyNo index)
        const int copyNo = hit->GetCopyNo();
        if (0 <= copyNo && copyNo < (int)fDE_MeV->size())
        {
          (*fDE_MeV)[copyNo] += hit->GetEdep() / MeV;
          // 「そのセグメントの代表pid」を入れる（今は hit の pid をそのまま）
          (*fPID)[copyNo] = hit->GetPDG();
        }
      }
    }

    // ---------------------
    // Target (ntupleId = 2)
    // ---------------------
    static G4int targetHCID = -1;
    if (targetHCID < 0)
    {
      targetHCID = G4SDManager::GetSDMpointer()->GetCollectionID("TargetHC");
    }

    auto *targetHC = static_cast<TargetHitsCollection *>(hce->GetHC(targetHCID));
    if (targetHC && targetHC->GetSize() > 0)
    {
      const auto *thit = (*targetHC)[0];
      const double eT = thit->GetEdep() / MeV;
      // (A) ターゲット内訳print（毎回出すのがダルい→フラグでON/OFF）
      if (fDoTargetTrace)
      {
        // もし TargetSD が「直近イベントの内訳」を持ってて出せるならここで呼ぶのが理想
        // 例: static_cast<TargetSD*>(...)->DumpLastEventTrace();
        // 今はあなたがTargetSD側でG4coutしてるなら、そっちも同様にフラグで抑制するのが良い
        G4cout << "[TARGET] eventID=" << eventID
               << " targetE=" << eT << " MeV\n";
      }

      // (B) ピーク窓だけ rndm 保存（ON/OFF切替）
      if (fSaveRndmOnPeak)
      {
        if (fPeakEmin_MeV < eT && eT < fPeakEmax_MeV)
        {
          std::ostringstream os;
          os << fRndmDir << "/" << fRndmPrefix << eventID << ".txt";
          G4Random::saveEngineStatus(os.str().c_str());

          G4cout << "[PEAK-CAND] eventID=" << eventID
                 << " targetE=" << eT << " MeV"
                 << " saved " << os.str().c_str() << G4endl;
        }
      }

      ana->FillNtupleIColumn(2, 0, eventID);
      ana->FillNtupleDColumn(2, 1, thit->GetTime() / ns);
      ana->FillNtupleDColumn(2, 2, thit->GetEdep() / MeV);
      ana->FillNtupleIColumn(2, 3, thit->GetPDG());
      ana->FillNtupleDColumn(2, 4, thit->GetMomentum().x() / MeV);
      ana->FillNtupleDColumn(2, 5, thit->GetMomentum().y() / MeV);
      ana->FillNtupleDColumn(2, 6, thit->GetMomentum().z() / MeV);
      ana->AddNtupleRow(2);
    }
  }

  // --------------------------
  // Event summary (ntupleId=0)
  // --------------------------
  ana->FillNtupleIColumn(0, 0, eventID);
  ana->FillNtupleDColumn(0, 1, fEvtEdepCell / MeV);
  ana->FillNtupleDColumn(0, 2, fEvtEdepTarget / MeV);
  ana->AddNtupleRow(0);

  // --------------------------
  // calarr (ntupleId=3)
  // --------------------------
  ana->FillNtupleIColumn(3, 0, eventID);
  ana->AddNtupleRow(3);
}