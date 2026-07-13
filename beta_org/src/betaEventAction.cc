#include "betaEventAction.hh"
#include "BetaConfig.hh"

#include "G4AnalysisManager.hh"
#include "G4Event.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include "G4SDManager.hh"
#include "G4GenericMessenger.hh"
#include "G4HCofThisEvent.hh"

#include "CalorimeterHit.hh"
#include "DetectorResponse.hh"
#include "HodoscopeHit.hh"
#include "TargetHit.hh"
#include "Constant.hh"

// #include "betaLineage.hh"

#include <algorithm>
#include <string>
#include <sstream>
#include <iostream>

betaEventAction::betaEventAction()
    : G4UserEventAction(),
      fDE_MeV(BetaConfig::Instance().NCells(), 0.0),
      fPID(BetaConfig::Instance().NCells(), 0),
      fTHDE_MeV(nSegTH, 0.0),
      fTHTime_ns(nSegTH, -1.0),
      fTHTimeLeft_ns(nSegTH, -1.0),
      fTHTimeRight_ns(nSegTH, -1.0),
      fTHTimeLeftMinusRight_ns(nSegTH, -9999.0),
      fTHZReco_mm(nSegTH, -9999.0),
      fTHPathTruth_mm(nSegTH, 0.0),
      fTLCDE_MeV(nSegTLC, 0.0),
      fTLCCherenkovTime_ns(nSegTLC, -1.0),
      fTLCPath_mm(nSegTLC, 0.0),
      fTLCCherenkovPath_mm(nSegTLC, 0.0),
      fTLCCherenkovExpectedPhotons(nSegTLC, 0.0)
{
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
  fEvtEdepPC = 0.;

  std::fill(fDE_MeV.begin(), fDE_MeV.end(), 0.0);
  std::fill(fPID.begin(), fPID.end(), 0);
  std::fill(fTHDE_MeV.begin(), fTHDE_MeV.end(), 0.0);
  std::fill(fTHTime_ns.begin(), fTHTime_ns.end(), -1.0);
  std::fill(fTHTimeLeft_ns.begin(), fTHTimeLeft_ns.end(), -1.0);
  std::fill(fTHTimeRight_ns.begin(), fTHTimeRight_ns.end(), -1.0);
  std::fill(fTHTimeLeftMinusRight_ns.begin(), fTHTimeLeftMinusRight_ns.end(), -9999.0);
  std::fill(fTHZReco_mm.begin(), fTHZReco_mm.end(), -9999.0);
  std::fill(fTHPathTruth_mm.begin(), fTHPathTruth_mm.end(), 0.0);
  std::fill(fTLCDE_MeV.begin(), fTLCDE_MeV.end(), 0.0);
  std::fill(fTLCCherenkovTime_ns.begin(), fTLCCherenkovTime_ns.end(), -1.0);
  std::fill(fTLCPath_mm.begin(), fTLCPath_mm.end(), 0.0);
  std::fill(fTLCCherenkovPath_mm.begin(), fTLCCherenkovPath_mm.end(), 0.0);
  std::fill(fTLCCherenkovExpectedPhotons.begin(), fTLCCherenkovExpectedPhotons.end(), 0.0);
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
    if (fCalHCID < 0)
      fCalHCID = G4SDManager::GetSDMpointer()->GetCollectionID("CalorimeterHC");

    auto *calHC = static_cast<CalorimeterHitsCollection *>(hce->GetHC(fCalHCID));
    if (calHC)
    {
      for (size_t i = 0; i < calHC->GetSize(); i++)
      {
        const auto *hit = (*calHC)[i];

        if (BetaConfig::Instance().WriteCalHit())
        {
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
        }

        // calarr vector (copyNo index)
        const int copyNo = hit->GetCopyNo();
        if (0 <= copyNo && copyNo < (int)fDE_MeV.size())
        {
          fDE_MeV[copyNo] += hit->GetEdep() / MeV;
          // 「そのセグメントの代表pid」を入れる（今は hit の pid をそのまま）
          fPID[copyNo] = hit->GetPDG();
        }
      }
    }

    // ---------------------
    // Target (ntupleId = 2)
    // ---------------------
    if (fTargetHCID < 0)
      fTargetHCID = G4SDManager::GetSDMpointer()->GetCollectionID("TargetHC");

    auto *targetHC = static_cast<TargetHitsCollection *>(hce->GetHC(fTargetHCID));
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

    if (fTHHCID < 0)
      fTHHCID = G4SDManager::GetSDMpointer()->GetCollectionID("THHC");
    auto *thHC = static_cast<HodoscopeHitsCollection *>(hce->GetHC(fTHHCID));
    if (thHC)
    {
      for (size_t i = 0; i < thHC->GetSize(); ++i)
      {
        const auto *hit = (*thHC)[i];
        const int copyNo = hit->GetCopyNo();
        if (copyNo < 0 || copyNo >= nSegTH)
          continue;
        fTHDE_MeV[copyNo] = hit->GetEdep() / MeV;
        if (hit->GetTime() < 1.0e29)
          fTHTime_ns[copyNo] = hit->GetTime() / ns;
        fTHPathTruth_mm[copyNo] = hit->GetChargedPath() / mm;

        const G4double leftTime = hit->GetTHLeftArrivalTime();
        const G4double rightTime = hit->GetTHRightArrivalTime();
        if (leftTime < 1.0e29 && rightTime < 1.0e29)
        {
          const G4double timeDifference = leftTime - rightTime;
          fTHTimeLeft_ns[copyNo] = leftTime / ns;
          fTHTimeRight_ns[copyNo] = rightTime / ns;
          fTHTimeLeftMinusRight_ns[copyNo] = timeDifference / ns;
          fTHZReco_mm[copyNo] =
              0.5 * DetectorResponse::THEffectiveLightSpeed * timeDifference / mm;
        }
      }
    }

    if (fTLCHCID < 0)
      fTLCHCID = G4SDManager::GetSDMpointer()->GetCollectionID("TLCHC");
    auto *tlcHC = static_cast<HodoscopeHitsCollection *>(hce->GetHC(fTLCHCID));
    if (tlcHC)
    {
      for (size_t i = 0; i < tlcHC->GetSize(); ++i)
      {
        const auto *hit = (*tlcHC)[i];
        const int copyNo = hit->GetCopyNo();
        if (copyNo < 0 || copyNo >= nSegTLC)
          continue;
        fTLCDE_MeV[copyNo] = hit->GetEdep() / MeV;
        if (hit->GetCherenkovTime() < 1.0e29)
          fTLCCherenkovTime_ns[copyNo] = hit->GetCherenkovTime() / ns;
        fTLCPath_mm[copyNo] = hit->GetChargedPath() / mm;
        fTLCCherenkovPath_mm[copyNo] = hit->GetCherenkovPath() / mm;
        fTLCCherenkovExpectedPhotons[copyNo] =
            hit->GetCherenkovExpectedPhotons();
      }
    }
  }

  // --------------------------
  // Event summary (ntupleId=0)
  // --------------------------
  ana->FillNtupleIColumn(0, 0, eventID);
  ana->FillNtupleDColumn(0, 1, fEvtEdepCell / MeV);
  ana->FillNtupleDColumn(0, 2, fEvtEdepTarget / MeV);
  ana->FillNtupleDColumn(0, 3, fEvtEdepPC / MeV);
  ana->AddNtupleRow(0);

  // --------------------------
  // calarr (ntupleId=3)
  // --------------------------
  ana->FillNtupleIColumn(3, 0, eventID);
  ana->AddNtupleRow(3);

  ana->FillNtupleIColumn(5, 0, eventID);
  ana->AddNtupleRow(5);
  ana->FillNtupleIColumn(6, 0, eventID);
  ana->AddNtupleRow(6);
}
