#ifndef CellParameterisation_H
#define CellParameterisation_H 1

#include "calConstant.hh"

#include "globals.hh"
#include "G4VPVParameterisation.hh"

#include <array>

class G4VPhysicalVolume;

/// Calorimeter cell parameterisation

class CellParameterisation : public G4VPVParameterisation
{
  public:
    CellParameterisation();
    virtual ~CellParameterisation();
    
    virtual void ComputeTransformation(
                   const G4int copyNo,G4VPhysicalVolume *physVol) const;
    
  private:
  //std::array<G4double, kNofEmCells> fXCell;
  //std::array<G4double, kNofEmCells> fYCell;
  std::array<G4double, 2> fXCell;
  std::array<G4double, 2> fYCell;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
