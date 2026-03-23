// calEventAction.hh
#ifndef calEventAction_h
#define calEventAction_h 1

//#include "ParamMan.hh"
#include "G4UserEventAction.hh"
#include "globals.hh"
#include "calConstant.hh"

class calEventAction : public G4UserEventAction
{
public:
  calEventAction();
  virtual ~calEventAction();
  //calEventAction(ParamManager*);

  virtual void  BeginOfEventAction(const G4Event* event);
  virtual void    EndOfEventAction(const G4Event* event);
    
  void AddAbs(G4double de, G4double dl);
    
private:
  //ParamManager* paramMan;

  G4int evID;
  G4double beamx,beamy,beamz;
  G4double beamMom;
  G4double beamE;
  G4double fEnergyAbs;
  G4double fTrackLAbs; 
};

// inline functions

inline void calEventAction::AddAbs(G4double de, G4double dl) {
  fEnergyAbs += de; 
  fTrackLAbs += dl;
}
                     
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
