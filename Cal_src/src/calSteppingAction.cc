// calSteppingAction.cc
#include "calSteppingAction.hh"
#include "calEventAction.hh"
#include "calDetectorConstruction.hh"

#include "G4Step.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

calSteppingAction::calSteppingAction(
                      const calDetectorConstruction* detectorConstruction,
                      calEventAction* eventAction)
  : G4UserSteppingAction(),
    fDetConstruction(detectorConstruction),
    fEventAction(eventAction)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

calSteppingAction::~calSteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void calSteppingAction::UserSteppingAction(const G4Step* step)
{
// Collect energy and track length step by step
  const G4StepPoint* preStepPoint = step->GetPreStepPoint();
  const G4StepPoint* postStepPoint = step->GetPostStepPoint();

  // get volume of the current step
  auto volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();

  // energy deposit
  auto edep = step->GetTotalEnergyDeposit();

  //G4String particleName = step->GetTrack()->GetDefinition()->GetParticleName();
  G4String proc = step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();

//  if(proc == "Decay"){
//    const std::vector<const G4Track*>* secondary = step->GetSecondaryInCurrentStep();
//    for(int lp=0; lp<(*secondary).size();lp++){
//      auto particle = (*secondary)[lp]->GetDefinition();
//      G4String pName = (*secondary)[lp]->GetDefinition()->GetParticleName();
//      if(pName=="proton"){ DecayFlag=1;}
//      if(pName=="neutron"){ DecayFlag=2;}
//      if(DecayFlag>0){ G4cout << "SSSSSSSSS"<< pName << G4endl;}
//    }
//    G4cout << "AAAA " << DecayFlag << G4endl;
//    getchar();
//  }

  // step length
  G4double stepLength = 0.;
  if ( step->GetTrack()->GetDefinition()->GetPDGCharge() != 0. ) {
    stepLength = step->GetStepLength();
  }
  if ( volume == fDetConstruction->GetCellPV() ) {
    //if ( volume == fDetConstruction->GetAbsorberPV() ) {
    //G4int Cellid = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber();
    fEventAction->AddAbs(edep,stepLength);
  }
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

