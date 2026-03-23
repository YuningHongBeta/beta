// calDetectorConstruction.cc
#include "calDetectorConstruction.hh"
#include "calConstant.hh"
#include "TargetSD.hh"
#include "CalorimeterSD.hh"
#include "ScintillatorSD.hh"
// #include "calBiasingOperator.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4MaterialTable.hh"
#include "G4NistManager.hh"

#include "G4Isotope.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal G4GlobalMagFieldMessenger *calDetectorConstruction::fMagFieldMessenger = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

calDetectorConstruction::calDetectorConstruction()
    : G4VUserDetectorConstruction(),
      fAbsorberPV(nullptr),
      fCellPV(nullptr),
      fVisAttributes(),
      fCheckOverlaps(true)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

calDetectorConstruction::~calDetectorConstruction()
{
  for (auto visAttributes : fVisAttributes)
  {
    delete visAttributes;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *calDetectorConstruction::Construct()
{
  // Define materials
  DefineMaterials();

  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void calDetectorConstruction::DefineMaterials()
{
  // Lead material defined using NIST Manager
  auto nistManager = G4NistManager::Instance();
  nistManager->FindOrBuildMaterial("G4_Pb");
  nistManager->FindOrBuildMaterial("G4_C");
  nistManager->FindOrBuildMaterial("G4_Li");
  nistManager->FindOrBuildMaterial("G4_BGO");
  nistManager->FindOrBuildMaterial("G4_CESIUM_IODIDE");
  nistManager->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
  nistManager->FindOrBuildMaterial("G4_AIR");
  nistManager->FindOrBuildMaterial("G4_PLEXIGLASS");

  G4double z;
  G4double a;
  G4double density;
  // Vacuum
  new G4Material("Galactic", z = 1., a = 1.01 * g / mole, density = universe_mean_density,
                 kStateGas, 2.73 * kelvin, 3.e-18 * pascal);
  // Liquid argon material
  new G4Material("liquidArgon", z = 18., a = 39.95 * g / mole, density = 1.390 * g / cm3);

  // Isotopes (atomic masses are in g/mole)
  auto Li6 = new G4Isotope("Li6", /*Z=*/3, /*N=*/6, 6.0151223 * g / mole);
  auto Li7 = new G4Isotope("Li7", /*Z=*/3, /*N=*/7, 7.0160040 * g / mole);

  // Element with enriched isotopic composition
  auto elLi_enr = new G4Element("Lithium_enriched", "Li_enr", 2);
  elLi_enr->AddIsotope(Li6, 90. * perCent); // atomic fraction
  elLi_enr->AddIsotope(Li7, 10. * perCent);

  // Material (single-element)
  auto rho = 0.47 * g / cm3; // <-- as requested
  auto matLi_enr = new G4Material("Li6_90pct", rho, 1, kStateSolid);
  matLi_enr->AddElement(elLi_enr, 1);

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *calDetectorConstruction::DefineVolumes()
{
  // Geometry parameters
  auto worldSizeXY = 3 * Rmax;
  auto worldSizeZ = 3 * Rmax;

  // Get materials
  auto defaultMaterial = G4Material::GetMaterial("Galactic");
  auto absorberMaterial = G4Material::GetMaterial("G4_CESIUM_IODIDE");
  auto ScintiMaterial = G4Material::GetMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
  auto TargetMaterial = G4Material::GetMaterial("G4_Li"); // coution! target should be Li6
  auto AirMaterial = G4Material::GetMaterial("G4_AIR");
  auto acrylMaterial = G4Material::GetMaterial("G4_PLEXIGLASS");
  auto LiTargetMaterial = G4Material::GetMaterial("Li6_90pct");

  //  if ( ! defaultMaterial || ! absorberMaterial ) {
  //    G4ExceptionDescription msg;
  //    msg << "Cannot retrieve materials already defined.";
  //    G4Exception("calDetectorConstruction::DefineVolumes()",
  //		"MyCode0001", FatalException, msg);
  //  }

  //
  // World
  //
  auto worldS = new G4Box("World",                                           // its name
                          worldSizeXY / 2, worldSizeXY / 2, worldSizeZ / 2); // its size

  auto worldLV = new G4LogicalVolume(
      worldS,          // its solid
      defaultMaterial, // its material
      "World");        // its name

  auto worldPV = new G4PVPlacement(
      0,               // no rotation
      G4ThreeVector(), // at (0,0,0)
      worldLV,         // its logical volume
      "World",         // its name
      0,               // its mother  volume
      false,           // no boolean operation
      0,               // copy number
      fCheckOverlaps); // checking overlaps

  //
  // Calorimeter //BOX
  //
  if (shapeFlag == 0)
  {
    auto calorimeterS = new G4Box("Calorimeter",                                        // its name
                                  calorSizeXY / 2, calorSizeXY / 2, absoThickness / 2); // its size

    calorLV = new G4LogicalVolume(
        calorimeterS,    // its solid
        defaultMaterial, // its material
        "Calorimeter");  // its name

    new G4PVPlacement(
        0,               // no rotation
        G4ThreeVector(), // at (0,0,0)
        calorLV,         // its logical volume
        "Calorimeter",   // its name
        worldLV,         // its mother  volume
        false,           // no boolean operation
        0,               // copy number
        fCheckOverlaps); // checking overlaps

    G4Box *cellS = new G4Box("Cell",                                                            // its name
                             calorSizeXY / 2, calorSizeXY / 2, absoThickness / (nb_cryst * 2)); // its size

    cellLV = new G4LogicalVolume(
        cellS,            // its solid
        absorberMaterial, // its material
        "Cell");          // its name

    G4int CopyNo = 0;
    for (G4int icrys = 0; icrys < nb_cryst; icrys++)
    {
      G4RotationMatrix rotm = G4RotationMatrix();
      G4double posZ = absoThickness * ((1. - nb_cryst + 2. * icrys) / (2. * nb_cryst));
      G4ThreeVector pos = G4ThreeVector(0, 0, posZ);
      G4Transform3D transform = G4Transform3D(rotm, pos);
      fCellPV = new G4PVPlacement(transform,       // rotation,position
                                  cellLV,          // its logical volume
                                  "cellPhysical",  // its name
                                  calorLV,         // its mother
                                  false,           // no boolean operation
                                  CopyNo,          // copy number
                                  fCheckOverlaps); // checking overlaps
      CopyNo++;
    } // for icrys
  }

  //
  // Calorimeter // Sphere
  //
  if (shapeFlag == 1)
  {
    G4Sphere *calorimeterS = new G4Sphere("Calorimeter",                // its name
                                          Rmin - ScintiThickness, Rmax, // Rmin, Rmax
                                          0 * deg, 360 * deg,           // phi
                                          0 * deg, 180 * deg);          // theta

    calorLV = new G4LogicalVolume(
        calorimeterS,    // its solid
        defaultMaterial, // its material
        "Calorimeter");  // its name

    new G4PVPlacement(
        0,               // no rotation
        G4ThreeVector(), // at (0,0,0)
        calorLV,         // its logical volume
        "Calorimeter",   // its name
        worldLV,         // its mother  volume
        false,           // no boolean operation
        0,               // copy number
        fCheckOverlaps); // checking overlaps

    // ==== ConstructCalorimeter ====: (innner radius, outer radius, material, cellLV, mother);
    G4String mat_Cal = "G4_BGO";
    G4String mat_Scinti = "G4_PLASTIC_SC_VINYLTOLUENE";
    // G4String mat_AC = """;

    // CsI(or BGO) crystals
    ConstructCalorimeter(Rmin, Rmax, mat_Cal, cellLV_Cal, calorLV);
    // plastic scincillator
    // ConstructCalorimeter(Rmin-ScintiThickness, Rmin, mat_Scinti, cellLV_Scinti, calorLV);

    // Target
    //  auto TargetS
    //    = new G4Box("Target",     // its name
    //    TargetSizeX/2, TargetSizeY/2, TargetThickness/2); // its size
    auto TargetS = new G4Tubs("Target",
                              0,
                              TargetRadius,
                              TargetLength / 2,
                              0 * deg,
                              360 * deg);

    TargetLV = new G4LogicalVolume(
        TargetS,          // its solid
        LiTargetMaterial, // its material
        "Target");        // its name
    new G4PVPlacement(
        0,               // no rotation
        G4ThreeVector(), // at (0,0,0)
        TargetLV,        // its logical volume
        "Target",        // its name
        worldLV,         // its mother  volume
        false,           // no boolean operation
        0,               // copy number
        fCheckOverlaps); // checking overlaps

    auto THmotherS = new G4Tubs("THmother",
                                ThRmin,
                                ThRmax,
                                TLength / 2,
                                0 * deg,
                                360 * deg);
    THmotherLV = new G4LogicalVolume(
        THmotherS,   // its solid
        AirMaterial, // its material
        "THmother"); // its name

    new G4PVPlacement(
        0,               // no rotation
        G4ThreeVector(), // at (0,0,0)
        THmotherLV,      // its logical volume
        "THmother",      // its name
        worldLV,         // its mother  volume
        false,           // no boolean operation
        0,               // copy number
        fCheckOverlaps); // checking overlaps

    auto THsegS = new G4Tubs("THseg",
                             ThRmin,
                             ThRmax,
                             TLength / 2,
                             0 * deg,
                             dPhiTH);
    THsegLV = new G4LogicalVolume(
        THsegS,         // its solid
        ScintiMaterial, // its material
        "THseg");       // its name
    new G4PVReplica(
        "THseg",    // its name
        THsegLV,    // its logical volume
        THmotherLV, // its mother volume
        kPhi,       // axis of replication
        nSegTH,     // number of replica
        dPhiTH);    // width of replica

    auto TLCmotherS = new G4Tubs("TLCmother",
                                 TlcRmin,
                                 TlcRmax,
                                 TLength / 2,
                                 0 * deg,
                                 360 * deg);
    TLCmotherLV = new G4LogicalVolume(
        TLCmotherS,   // its solid
        AirMaterial,  // its material
        "TLCmother"); // its name
    new G4PVPlacement(
        0,               // no rotation
        G4ThreeVector(), // at (0,0,0)
        TLCmotherLV,     // its logical volume
        "TLCmother",     // its name
        worldLV,         // its mother  volume
        false,           // no boolean operation
        0,               // copy number
        fCheckOverlaps); // checking overlaps
    auto TLCsegS = new G4Tubs("TLCseg",
                              TlcRmin,
                              TlcRmax,
                              TLength / 2,
                              0 * deg,
                              dPhiTLC);
    TLCsegLV = new G4LogicalVolume(
        TLCsegS,       // its solid
        acrylMaterial, // its material
        "TLCseg");     // its name
    new G4PVReplica(
        "TLCseg",    // its name
        TLCsegLV,    // its logical volume
        TLCmotherLV, // its mother volume
        kPhi,        // axis of replication
        nSegTLC,     // number of replica
        dPhiTLC);    // width of replica
  }

  //
  // print parameters
  //
  //  G4cout
  //    << G4endl
  //    << "------------------------------------------------------------" << G4endl
  //    << "---> The calorimeter is " << nofLayers << " layers of: [ "
  //    << absoThickness/mm << "mm of " << absorberMaterial->GetName()
  //    << " + "
  //    << gapThickness/mm << "mm of " << gapMaterial->GetName() << " ] " << G4endl
  //    << "------------------------------------------------------------" << G4endl;

  //
  // Visualization attributes
  //
  worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());
  THmotherLV->SetVisAttributes(G4VisAttributes::GetInvisible());
  TLCmotherLV->SetVisAttributes(G4VisAttributes::GetInvisible());

  auto visAttributes = new G4VisAttributes(G4Colour::Green());
  visAttributes->SetVisibility(false);
  calorLV->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);

  visAttributes = new G4VisAttributes(G4Colour::Blue());
  visAttributes->SetVisibility(false);
  worldLV->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);

  // Calorimeter //BOX
  if (shapeFlag == 0)
  {
    visAttributes = new G4VisAttributes(G4Colour::Cyan());
    visAttributes->SetVisibility(true);
    cellLV->SetVisAttributes(visAttributes);
  }

  //  // Calorimeter // Sphere
  if (shapeFlag == 1)
  {
    visAttributes = new G4VisAttributes(G4Colour::Cyan());
    visAttributes->SetVisibility(true);
    for (G4int j = 0; j < nLayer; j++)
    {
      cellLV_Cal[j]->SetVisAttributes(visAttributes);
    }
    fVisAttributes.push_back(visAttributes);

    // visAttributes = new G4VisAttributes(G4Colour::Blue());
    // visAttributes->SetVisibility(true);
    // for (G4int j = 0; j < nLayer  ; j++){
    //   cellLV_Scinti[j]->SetVisAttributes(visAttributes);
    // }
    // fVisAttributes.push_back(visAttributes);

    visAttributes = new G4VisAttributes(G4Colour::White());
    visAttributes->SetVisibility(true);
    visAttributes->SetForceSolid(true);
    TargetLV->SetVisAttributes(visAttributes);
    fVisAttributes.push_back(visAttributes);

    visAttributes = new G4VisAttributes(G4Colour::Red());
    visAttributes->SetVisibility(true);
    THsegLV->SetVisAttributes(visAttributes);
    fVisAttributes.push_back(visAttributes);

    visAttributes = new G4VisAttributes(G4Colour::Green());
    visAttributes->SetVisibility(true);
    TLCsegLV->SetVisAttributes(visAttributes);
    fVisAttributes.push_back(visAttributes);
  }

  //
  // Always return the physical World
  //
  return worldPV;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LogicalVolume *calDetectorConstruction::ConstructCell(G4double rmin, G4double rmax, G4double center, G4String mat)
{
  auto Material = G4Material::GetMaterial(mat);

  G4Sphere *cellS = new G4Sphere("Cell",                                                  // its name
                                 rmin, rmax,                                              // Rmin, Rmax
                                 0 * deg, phi * deg,                                      // phi
                                 center * deg - (cellAngle * deg / 2.), cellAngle * deg); // theta

  cellLV = new G4LogicalVolume(
      cellS,    // its solid
      Material, // its material
      "Cell");  // its name

  return cellLV;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Transform3D calDetectorConstruction::Rotate(G4double icrys)
{
  G4ThreeVector pos = G4ThreeVector(0, 0, 0);
  G4RotationMatrix rotm = G4RotationMatrix();
  rotm.rotateZ(phi * icrys * deg);
  G4Transform3D transform = G4Transform3D(rotm, pos);
  return transform;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void calDetectorConstruction::ConstructCalorimeter(G4double rmin, G4double rmax, G4String mat, G4LogicalVolume *CELL[], G4LogicalVolume *mother)
{
  G4int CopyNo = 0;
  G4double center = 90;

  for (G4int j = 0; j < nLayer; j++)
  {
    // center = (180 - (1 - nLayer + 2 * j) * cellAngle) / 2.;
    // CELL[j] = ConstructCell(rmin, rmax, center, mat);
    center = thetaMin + (j + 0.5) * cellAngle; // [deg]
    CELL[j] = ConstructCell(rmin, rmax, center, mat);
    for (G4int icrys = 0; icrys < nSector; icrys++)
    {
      G4Transform3D transform = Rotate(icrys);
      fCellPV = new G4PVPlacement(transform, CELL[j], "cellPhysical",
                                  mother, false, CopyNo, fCheckOverlaps);
      CopyNo++;
    }
  } // for j
}

void calDetectorConstruction::ConstructSDandField()
{

  // sensitive detectors -----------------------------------------------------
  auto sdManager = G4SDManager::GetSDMpointer();
  G4String SDname;

  auto calorimeter = new CalorimeterSD(SDname = "/calorimeter");
  sdManager->AddNewDetector(calorimeter);
  if (shapeFlag == 0)
  {
    cellLV->SetSensitiveDetector(calorimeter);
    // calBiasingOperator *biasOp = new calBiasingOperator("CellBiasOp", 1.0);
    // biasOp->AttachTo(cellLV);
  }
  if (shapeFlag == 1)
  {
    auto target = new TargetSD(SDname = "/target");
    sdManager->AddNewDetector(target);
    TargetLV->SetSensitiveDetector(target);

    for (G4int j = 0; j < nLayer; j++)
    {
      cellLV_Cal[j]->SetSensitiveDetector(calorimeter);
    }

    // auto scintillator = new ScintillatorSD(SDname="/scintillator");
    // sdManager->AddNewDetector(scintillator);
    // for (G4int j = 0; j < nLayer  ; j++){
    //   cellLV_Scinti[j]->SetSensitiveDetector(scintillator);
    // }
  }

  // magnetic field ----------------------------------------------------------
  //  G4ThreeVector fieldValue;
  //  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  //  fMagFieldMessenger->SetVerboseLevel(1);
  //  Register the field messenger for deleting
  //  G4AutoDelete::Register(fMagFieldMessenger);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
