// TargetSD.cc
#include "TargetSD.hh"
#include "calEventAction.hh"
#include "calConstant.hh"

#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

#include "TROOT.h"
#include "TH1.h"
#include "TTree.h"
#include "TFile.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TargetSD::TargetSD(G4String name)
  : G4VSensitiveDetector(name)
{
  collectionName.insert("TargetColl");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TargetSD::~TargetSD()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TargetSD::Initialize(G4HCofThisEvent* hce)
{
  etot_Target=0.;
  edepbuf_Target =0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool TargetSD::ProcessHits(G4Step*step, G4TouchableHistory*)
{
  const G4StepPoint* preStepPoint = step->GetPreStepPoint();
  const G4StepPoint* postStepPoint = step->GetPostStepPoint();
 
  auto edep = step->GetTotalEnergyDeposit();
  if (edep==0.) return true;
  
  auto touchable = step->GetPreStepPoint()->GetTouchable();
  auto detector = touchable->GetVolume();

  // get detector Name 
  G4String detectorName = detector->GetName();    
  G4String particleName = step->GetTrack()->GetDefinition()->GetParticleName();

  edepbuf_Target += edep;  

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TargetSD::EndOfEvent(G4HCofThisEvent* hce)
{
  etot_Target = edepbuf_Target;
  
  TTree* tree = (TTree*)gROOT-> FindObject("tree");
  tree->SetBranchAddress("dEtot_Target"	, &etot_Target);

}

