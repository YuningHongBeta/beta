// TargetSD.hh
#ifndef TargetSD_h
#define TargetSD_h 1

#include "G4VSensitiveDetector.hh"
#include "calConstant.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class TargetSD : public G4VSensitiveDetector
{   
  public:
    TargetSD(G4String name);
    virtual ~TargetSD();
    
    virtual void Initialize(G4HCofThisEvent*HCE);
    virtual void EndOfEvent(G4HCofThisEvent*HCE);
    virtual G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
    
  private:
    G4double etot_Target;
    G4double edepbuf_Target;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
