// calActionInitialization.cc

#include "calActionInitialization.hh"
#include "calPrimaryGeneratorAction.hh"
#include "calRunAction.hh"
#include "calEventAction.hh"
#include "calSteppingAction.hh"
#include "calDetectorConstruction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

calActionInitialization::calActionInitialization
                            (calDetectorConstruction* detConstruction)
 : G4VUserActionInitialization(),
   fDetConstruction(detConstruction)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

calActionInitialization::~calActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void calActionInitialization::BuildForMaster() const
{
  SetUserAction(new calRunAction);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void calActionInitialization::Build() const
{
  SetUserAction(new calPrimaryGeneratorAction);
  SetUserAction(new calRunAction);
  auto eventAction = new calEventAction;
  SetUserAction(eventAction);
  SetUserAction(new calSteppingAction(fDetConstruction,eventAction));
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
