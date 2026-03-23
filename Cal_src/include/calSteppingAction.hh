// calSteppingAction.hh
#ifndef calSteppingAction_h
#define calSteppingAction_h 1

#include "G4UserSteppingAction.hh"

class calDetectorConstruction;
class calEventAction;

class calSteppingAction : public G4UserSteppingAction
{
public:
  calSteppingAction(const calDetectorConstruction* detectorConstruction,
                    calEventAction* eventAction);
  virtual ~calSteppingAction();

  virtual void UserSteppingAction(const G4Step* step);
    
private:
  const calDetectorConstruction* fDetConstruction;
  calEventAction*  fEventAction;  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
