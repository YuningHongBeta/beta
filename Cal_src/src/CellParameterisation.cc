#include "CellParameterisation.hh"
#include "calConstant.hh"

#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CellParameterisation::CellParameterisation()
: G4VPVParameterisation()
{
  for (auto copyNo=0; copyNo<nb_cryst; copyNo++) {
    auto column = 1;//copyNo / kNofEmRows ;
    auto row = 1;//copyNo % kNofEmRows;
    fXCell[copyNo] = (column-9)*15.*cm - 7.5*cm;
    fYCell[copyNo] = (row-1)*15*cm - 7.5*cm;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CellParameterisation::~CellParameterisation()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CellParameterisation::ComputeTransformation(
       const G4int copyNo,G4VPhysicalVolume *physVol) const
{
  physVol->SetTranslation(G4ThreeVector(fXCell[copyNo],fYCell[copyNo],0.));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
