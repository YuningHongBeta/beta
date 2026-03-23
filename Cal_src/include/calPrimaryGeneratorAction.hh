/// \file calPrimaryGeneratorAction.hh
/// \brief Definition of the calPrimaryGeneratorAction class

#ifndef calPrimaryGeneratorAction_h
#define calPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
//#include "ParamMan.hh"

class G4ParticleGun;
class G4Event;


class calPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  calPrimaryGeneratorAction();    
  virtual ~calPrimaryGeneratorAction();
  //calPrimaryGeneratorAction(ParamManager*);
  
  virtual void GeneratePrimaries(G4Event* event);

  virtual void init();
  virtual G4double incident(G4double ma, G4double pa, G4double ua, G4double va,
			    G4double mb, G4double pb, G4double ub, G4double vb);
  virtual G4double ThreeBodyDecay(G4double m1, G4double m2, G4double m3, G4double Ecm);

  // set methods
  void SetRandomFlag(G4bool value);

private:
  G4ParticleGun*  fParticleGun; // G4 particle gun
  //ParamManager* paramMan;
  G4double Beamx,Beamy,Beamz;
  G4double BeamPx,BeamPy,BeamPz;
  G4double BeamMom;
  G4double BeamE;
  G4int    decayFlag;

  G4double ma, pa, ua, va;
  G4double pa_x,pa_y,pa_z;
  G4double mb, pb, ub, vb;
  G4double pb_x,pb_y,pb_z;
  G4ThreeVector pA, pB, psum;
  G4double Ea, Eb, psum1, Esum1, Ecm; 
  G4ThreeVector beta;
  G4double E1cm, E2cm;
  G4double p1cm, p3;
  G4double m1, m2, m3, M;
  G4double rnd_Zm3, rnd_phim3;
  G4double pm3x,pm3y,pm3z;
  G4ThreeVector p12cm, p3cm;
  G4double m1m2, Mm3, m12, m23, m23min, m23max;
  G4double E12cm, E3cm, E2, E3, tmp1, tmp2, tmp3, Ee;
  
  G4ThreeVector beamPos;
  G4ThreeVector direction;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
