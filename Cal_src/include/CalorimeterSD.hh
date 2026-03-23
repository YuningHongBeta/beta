// CalorimeterSD.hh
#ifndef CalorimeterSD_h
#define CalorimeterSD_h 1

#include "G4VSensitiveDetector.hh"

#include "CalorimeterHit.hh"
#include "calConstant.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class CalorimeterSD : public G4VSensitiveDetector
{   
  public:
    CalorimeterSD(G4String name);
    virtual ~CalorimeterSD();
    
    virtual void Initialize(G4HCofThisEvent*HCE);
    virtual void EndOfEvent(G4HCofThisEvent*HCE);
    virtual G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
    
  private:
    CalorimeterHitsCollection* fHitsCollection;
    G4int fHCID;
    G4int cellNoFlag[nb_cryst];
    G4int multiHit;
    G4double etot;
    G4double edepbuf[nb_cryst];
    G4double firsthitposXbuf[nb_cryst];
    G4double firsthitposYbuf[nb_cryst];
    G4double firsthitposZbuf[nb_cryst];
    G4double firsthitTimebuf[nb_cryst];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
