#include "TargetSD.hh"

#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4StepPoint.hh"
#include "G4TouchableHistory.hh"
#include "G4VProcess.hh"

#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
#include "G4ios.hh"

#include <algorithm>
#include <vector>

TargetSD::TargetSD(const G4String& name, const G4String& hitsCollectionName)
  : G4VSensitiveDetector(name)
{
  collectionName.insert(hitsCollectionName);
}

void TargetSD::Initialize(G4HCofThisEvent* hce)
{
  // --- hits collection ---
  fHitsCollection = new TargetHitsCollection(SensitiveDetectorName, collectionName[0]);

  if (fHCID < 0) {
    fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection);
  }
  hce->AddHitsCollection(fHCID, fHitsCollection);

  // --- reset per-event accumulators ---
  fEdepTot = 0.0;
  fEdepByPDG.clear();
  fEdepByProc.clear();

  // --- create one hit per event (index 0) ---
  fHit = new TargetHit();
  // timeは「そのvolumeで最初にエネルギーが落ちた時刻」を取りたいので初期は大きく
  fHit->SetTime(1e30);
  fHit->SetEdep(0.0);
  fHit->SetPDG(0);
  fHit->SetMomentum(G4ThreeVector(0,0,0));

  fHitsCollection->insert(fHit);
}

G4bool TargetSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
  const G4double edep = step->GetTotalEnergyDeposit();
  if (edep <= 0.0) return false;

  const G4StepPoint* pre  = step->GetPreStepPoint();
  const G4StepPoint* post = step->GetPostStepPoint();
  const G4Track* tr       = step->GetTrack();

  const G4double time = pre->GetGlobalTime();
  const G4int pdg     = tr->GetDefinition()->GetPDGEncoding();
  const G4ThreeVector mom = pre->GetMomentum();

  // --- process name ---
  const G4VProcess* p = nullptr;
  if (post) p = post->GetProcessDefinedStep();
  const std::string proc = (p ? p->GetProcessName() : "Undefined");

  // --- accumulate for tracing ---
  fEdepTot += edep;
  fEdepByPDG[pdg] += edep;
  fEdepByProc[proc] += edep;

  // --- update 1 hit (targetHC[0]) ---
  if (fHit) {
    // total Edep
    fHit->SetEdep(fHit->GetEdep() + edep);

    // earliest time wins
    if (time < fHit->GetTime()) {
      fHit->SetTime(time);
      fHit->SetPDG(pdg);
      fHit->SetMomentum(mom);
    }
  }

  return true;
}

void TargetSD::EndOfEvent(G4HCofThisEvent*)
{
  const G4double eT_MeV = fEdepTot / MeV;
  if (!(fPeakLoMeV < eT_MeV && eT_MeV < fPeakHiMeV)) return;

  // eventID も一緒に出したい（EndOfEventにはevtが来ないので RunManager から取る）
  G4int eventID = -1;
  if (auto* evt = G4RunManager::GetRunManager()->GetCurrentEvent()) {
    eventID = evt->GetEventID();
  }

  // G4cout << "\n=== [TARGET-TRACE] eventID=" << eventID
  //        << " targetE=" << eT_MeV << " MeV ===\n";

  // ---- by PDG (descending) ----
  {
    std::vector<std::pair<G4int, G4double>> v;
    v.reserve(fEdepByPDG.size());
    for (const auto& kv : fEdepByPDG) v.push_back(kv);

    std::sort(v.begin(), v.end(),
              [](const auto& a, const auto& b){ return a.second > b.second; });

  //   G4cout << "  -- by PDG (MeV) --\n";
  //   for (const auto& it : v) {
  //     G4cout << "    PDG " << it.first << " : " << (it.second/MeV) << " MeV\n";
  //   }
  }

  // ---- by process (descending) ----
  {
    std::vector<std::pair<std::string, G4double>> v;
    v.reserve(fEdepByProc.size());
    for (const auto& kv : fEdepByProc) v.push_back(kv);

    std::sort(v.begin(), v.end(),
              [](const auto& a, const auto& b){ return a.second > b.second; });

    // G4cout << "  -- by process (MeV) --\n";
    // for (const auto& it : v) {
    //   G4cout << "    " << it.first << " : " << (it.second/MeV) << " MeV\n";
    // }
  }

  // G4cout << "=== [TARGET-TRACE END] ===\n\n";
}