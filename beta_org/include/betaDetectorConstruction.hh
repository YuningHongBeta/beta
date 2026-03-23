// betaDetectorConstruction.hh
#ifndef betaDetectorConstruction_h
#define betaDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4Transform3D.hh"
#include "globals.hh"
#include "Constant.hh"

class G4VPhysicalVolume;
class G4VSensitiveDetector;
class G4GlobalMagFieldMessenger;
class G4VisAttributes;

class betaDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    betaDetectorConstruction();
    virtual ~betaDetectorConstruction();

  public:
    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();

    // get methods
    //
    const G4VPhysicalVolume* GetAbsorberPV() const;
    const G4VPhysicalVolume* GetCellPV() const;
  
  private:
    // methods
    //
    void DefineMaterials();
    G4VPhysicalVolume* DefineVolumes();
    void ConstructCalorimeter(G4double rmin, G4double rmax, G4String mat, G4LogicalVolume* CELL[], G4LogicalVolume* mother);
    G4LogicalVolume* ConstructCell(G4double rmin, G4double rmax, G4double center, G4String mat);
    G4Transform3D Rotate(G4double icrys);
  
    // data members
    //

    static G4ThreadLocal G4GlobalMagFieldMessenger*  fMagFieldMessenger; 
    G4LogicalVolume* calorLV;
    G4LogicalVolume* cellLV;
    G4LogicalVolume* cellLV_Cal[nLayer];     
    G4LogicalVolume* cellLV_Scinti[nLayer];     
    G4LogicalVolume* cellLV_AC[nLayer];     
    G4LogicalVolume* TargetLV;
    G4LogicalVolume* THsegLV;
    G4LogicalVolume* TLCsegLV;
    G4LogicalVolume* THmotherLV;
    G4LogicalVolume* TLCmotherLV;

    G4VPhysicalVolume*   fAbsorberPV; // the absorber physical volume
    G4VPhysicalVolume*   fCellPV;     // the Cell physical volume
    G4VPhysicalVolume*   fScintiPV;     // the Scintillator physical volume
    G4VPhysicalVolume*   fTargetPV;     // the Target physical volume

    
    std::vector<G4VisAttributes*> fVisAttributes;

    G4bool  fCheckOverlaps; // option to activate checking of volumes overlaps
};

// inline functions

inline const G4VPhysicalVolume* betaDetectorConstruction::GetAbsorberPV() const { 
  return fAbsorberPV; 
}
inline const G4VPhysicalVolume* betaDetectorConstruction::GetCellPV() const { 
  return fCellPV; 
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

