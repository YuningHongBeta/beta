// src/CalorimeterSD.cc  (MT-safe, quiet)

#include "CalorimeterSD.hh"
#include "CalorimeterHit.hh"
#include "betaStackingAction.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHandle.hh"
#include "G4StepPoint.hh"
#include "G4VProcess.hh"

CalorimeterSD::CalorimeterSD(const G4String& name,
                             const G4String& hitsCollectionName)
: G4VSensitiveDetector(name)
{
  collectionName.insert(hitsCollectionName);
}

void CalorimeterSD::Initialize(G4HCofThisEvent* hce)
{
  fHitsCollection =
    new CalorimeterHitsCollection(SensitiveDetectorName, collectionName[0]);

  if (fHCID < 0) {
    fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection);
  }
  hce->AddHitsCollection(fHCID, fHitsCollection);

  // event start reset
  fHitByCopyNo.clear();
}

G4bool CalorimeterSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
  if (!step) return false;

  const auto edep = step->GetTotalEnergyDeposit();
  if (edep <= 0.) return false;

  const auto* pre   = step->GetPreStepPoint();
  const auto* track = step->GetTrack();
  if (!pre || !track) return false;

  const auto th = pre->GetTouchableHandle();
  if (!th) return false;

  const int copyNo = th->GetCopyNumber();
  if (copyNo < 0) return false;

  const auto* cp = track->GetCreatorProcess();
  const G4String creatorName = cp ? cp->GetProcessName() : "Primary";

  // ensure vector size
  if ((size_t)copyNo >= fHitByCopyNo.size()) {
    fHitByCopyNo.resize(copyNo + 1, nullptr);
  }

  const auto time = pre->GetGlobalTime();
  const auto pdg  = track->GetDefinition()->GetPDGEncoding(); // first-hit pid
  const auto mom  = pre->GetMomentum();

  auto*& hit = fHitByCopyNo[copyNo];

  if (!hit) {
    // ---- first hit in this segment ----
    hit = new CalorimeterHit();
    hit->SetCopyNo(copyNo);

    hit->SetTime(time);     // first time
    hit->SetEdep(edep);     // start sum
    hit->SetPDG(pdg);       // first-hit PDG (fixed)
    hit->SetMomentum(mom);  // first momentum (fixed)
    hit->SetCreator(creatorName); // creator process name
    
    // Set Origin Type from Track Info
    if (auto* info = dynamic_cast<betaTrackInfo*>(track->GetUserInformation())) {
        hit->SetOriginType(info->GetOrigin());
    } else {
        hit->SetOriginType(0);
    }

    fHitsCollection->insert(hit);
  } else {
    // ---- accumulate only ----
    hit->SetEdep(hit->GetEdep() + edep);
  }

  return true;
}