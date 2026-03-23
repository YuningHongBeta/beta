/*                                                                                                                                  
  ParamManger.hh                                                                                                                    
*/

#ifndef ParamManager_h
#define ParamManager_h 1

//#include<TFile.h>                                                                                                                 
#include "G4String.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"
#include "calConstant.hh"
#include <fstream>
using namespace std;



class ParamManager
{
public:
  ParamManager();
  virtual ~ParamManager();

public:
  G4bool readparam();

private:
 // G4double beamx , beamy , beamz;
 // G4double beampx , beampy , beampz;
 // G4double beamMom;
 // G4double Ee;


public:
  G4double particle_energy = -100.;

  void test(G4double aEnergy){
    //particle_energy = aEnergy;
    G4cout <<  "EEEEEEe  "<< particle_energy << G4endl;
    getchar();
  }

 // void SetBeamPos(G4ThreeVector bpos){
 //   beamx = bpos.x();
 //   beamy = bpos.y();
 //   beamz = bpos.z();
 // }
 // void SetBeamMomdir(G4ThreeVector bmom){
 //   beampx = bmom.x();
 //   beampy = bmom.y();
 //   beampz = bmom.z();
 // }
 // void SetBeamMom(G4double mom){
 //   beamMom = mom;
 // }
 // void SetBeamE(G4double energy){
 //   //BeamE = energy;
 // }

public:
 // G4double GetBeamPosx() { return beamx; }
 // G4double GetBeamPosy() { return beamy; }
 // G4double GetBeamPosz() { return beamz; }
 // G4double GetBeamMomx() { return beampx; }
 // G4double GetBeamMomy() { return beampy; }
 // G4double GetBeamMomz() { return beampz; }
 // G4double GetBeamMom() { return beamMom; }
  G4double GetTest() { return particle_energy; }


};

#endif
