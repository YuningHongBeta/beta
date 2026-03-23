// ScintillatorSD.hh
#ifndef ScintillatorSD_h
#define ScintillatorSD_h 1

#include "G4VSensitiveDetector.hh"
#include "calConstant.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class ScintillatorSD : public G4VSensitiveDetector
{   
  public:
    ScintillatorSD(G4String name);
    virtual ~ScintillatorSD();
    
    virtual void Initialize(G4HCofThisEvent*HCE);
    virtual void EndOfEvent(G4HCofThisEvent*HCE);
    virtual G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
    
  private:
    //ScintillatorHitsCollection* fHitsCollection;
    //G4int fHCID;
    G4double etot_Scinti;
    G4double edepbuf_Scinti[nb_cryst];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
