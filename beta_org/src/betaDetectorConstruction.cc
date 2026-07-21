// betaDetectorConstruction.cc
#include "betaDetectorConstruction.hh"
#include "BGOeggGeometry.hh"
#include "BetaConfig.hh"
#include "Constant.hh"
#include "TargetSD.hh"
#include "CalorimeterSD.hh"
#include "HodoscopeSD.hh"

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
#include "G4Cons.hh"
#include "G4GenericTrap.hh"
#include "G4Sphere.hh"
#include "G4TwoVector.hh"
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

#include "betaBiasingOperator.hh"

#include <cmath>
#include <stdexcept>
#include <string>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal G4GlobalMagFieldMessenger *betaDetectorConstruction::fMagFieldMessenger = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

betaDetectorConstruction::betaDetectorConstruction()
    : G4VUserDetectorConstruction(),
      calorLV(nullptr),
      fAbsorberPV(nullptr),
      fCellPV(nullptr),
      fVisAttributes(),
      fCheckOverlaps(true)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

betaDetectorConstruction::~betaDetectorConstruction()
{
  for (auto visAttributes : fVisAttributes)
  {
    delete visAttributes;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *betaDetectorConstruction::Construct()
{
  // Define materials
  DefineMaterials();

  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void betaDetectorConstruction::DefineMaterials()
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
  auto C13 = new G4Isotope("C13", /*Z=*/6, /*N=*/13, 13.003354835 * g / mole);

  // Element with enriched isotopic composition
  auto elLi_enr = new G4Element("Lithium_enriched", "Li_enr", 2);
  elLi_enr->AddIsotope(Li6, 90. * perCent); // atomic fraction
  elLi_enr->AddIsotope(Li7, 10. * perCent);

  // Material (single-element)
  const auto &config = BetaConfig::Instance();
  auto rho = config.TargetDensityGCM3() * g / cm3;
  auto matLi_enr = new G4Material("Li6_90pct", rho, 1, kStateSolid);
  matLi_enr->AddElement(elLi_enr, 1);

  auto elC13 = new G4Element("Carbon13_enriched", "C13_enr", 1);
  elC13->AddIsotope(C13, 100. * perCent);
  auto matC13 = new G4Material("C13_100pct", rho, 1, kStateSolid);
  matC13->AddElement(elC13, 1);

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *betaDetectorConstruction::DefineVolumes()
{
  const auto &config = BetaConfig::Instance();
  const G4int nLayer = config.NLayer();
  const G4int nSector = config.NSector();
  const G4int nCells = config.NCells();
  const G4double calorRMin = config.RMinCm() * cm;
  const G4double calorThickness = config.ThicknessCm() * cm;
  const G4double calorRMax = calorRMin + calorThickness;

  // Geometry parameters
  auto worldSizeXY = config.BGOeggFrustum() ? 180.0 * cm : 3 * calorRMax;
  auto worldSizeZ = config.BGOeggFrustum() ? 180.0 * cm : 3 * calorRMax;

  // Get materials
  auto defaultMaterial = G4Material::GetMaterial("Galactic");
  auto absorberMaterial = G4Material::GetMaterial("G4_CESIUM_IODIDE");
  auto ScintiMaterial = G4Material::GetMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
  // auto TargetMaterial = G4Material::GetMaterial("G4_Li"); // coution! target should be Li6
  auto AirMaterial = G4Material::GetMaterial("G4_AIR");
  auto acrylMaterial = G4Material::GetMaterial("G4_PLEXIGLASS");
  auto TargetMaterial = G4Material::GetMaterial(config.TargetMaterial());

  //  if ( ! defaultMaterial || ! absorberMaterial ) {
  //    G4ExceptionDescription msg;
  //    msg << "Cannot retrieve materials already defined.";
  //    G4Exception("betaDetectorConstruction::DefineVolumes()",
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
                             calorSizeXY / 2, calorSizeXY / 2, absoThickness / (nCells * 2)); // its size

    cellLV = new G4LogicalVolume(
        cellS,            // its solid
        absorberMaterial, // its material
        "Cell");          // its name

    G4int CopyNo = 0;
    for (G4int icrys = 0; icrys < nCells; icrys++)
    {
      G4RotationMatrix rotm = G4RotationMatrix();
      G4double posZ = absoThickness * ((1. - nCells + 2. * icrys) / (2. * nCells));
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
    G4String mat_Cal = "G4_BGO";
    G4String mat_Scinti = "G4_PLASTIC_SC_VINYLTOLUENE";
    // G4String mat_AC = """;

    if (config.BGOeggFrustum())
    {
      // Published BGOegg frusta can extend beyond R=420 mm.  Place them
      // directly in the world so a vacuum mother does not overlap the target,
      // TH, or TLC volumes.
      ConstructBGOeggFrusta(
          mat_Cal, cellLV_Cal, worldLV, config.BgoZOffsetCm() * cm);
    }
    else
    {
      G4Sphere *calorimeterS = new G4Sphere(
          "Calorimeter", calorRMin - ScintiThickness, calorRMax,
          0 * deg, 360 * deg, 0 * deg, 180 * deg);

      calorLV = new G4LogicalVolume(
          calorimeterS, defaultMaterial, "Calorimeter");

      new G4PVPlacement(
          0, G4ThreeVector(0., 0., config.BgoZOffsetCm() * cm),
          calorLV, "Calorimeter", worldLV, false, 0, fCheckOverlaps);

      ConstructCalorimeter(
          calorRMin, calorRMax, mat_Cal, cellLV_Cal, calorLV);
    }
    // plastic scincillator
    // ConstructCalorimeter(Rmin-ScintiThickness, Rmin, mat_Scinti, cellLV_Scinti, calorLV);

    // Target
    //  auto TargetS
    //    = new G4Box("Target",     // its name
    //    TargetSizeX/2, TargetSizeY/2, TargetThickness/2); // its size
    auto TargetS = new G4Tubs("Target",
                              0,
                              config.TargetRadiusMm() * mm,
                              config.TargetLengthMm() * mm / 2,
                              0 * deg,
                              360 * deg);

    TargetLV = new G4LogicalVolume(
        TargetS,          // its solid
        TargetMaterial,   // its material
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

    if (config.PhotonCounterMode() != 0)
    {
      // Pb/plastic endcap annuli.  For the exact BGOegg they are required by
      // BetaConfig to start beyond the outermost frustum in z; their angular
      // range is held constant through the stack with conical boundaries.
      const G4int pcNLayers = config.PCNLayers();
      const G4double pcPbThickness = config.PCPbThicknessMm() * mm;
      const G4double pcScintiThickness = config.PCScintiThicknessMm() * mm;
      const G4double pcZFront = config.PCZFrontCm() * cm;
      const G4double pcLayerThickness = pcPbThickness + pcScintiThickness;
      const G4double pcZBack = pcZFront + pcNLayers * pcLayerThickness;
      auto pcPbVis = new G4VisAttributes(G4Colour::Gray());
      pcPbVis->SetVisibility(true);
      fVisAttributes.push_back(pcPbVis);
      auto pcScintiVis = new G4VisAttributes(G4Colour::Yellow());
      pcScintiVis->SetVisibility(true);
      fVisAttributes.push_back(pcScintiVis);

      // Legacy PCDown means the +z endcap.  The incident beam travels +z ->
      // -z, so this is physically upstream (the small BGOC opening).
      if (config.HasDownstreamPhotonCounter())
      {
        const G4double tanInner =
            std::tan(config.PCDownThetaInnerDeg() * deg);
        const G4double tanOuter =
            std::tan(config.PCDownThetaOuterDeg() * deg);
        for (G4int layer = 0; layer < pcNLayers; ++layer)
        {
          const G4double leadFront = pcZFront + layer * pcLayerThickness;
          const G4double leadBack = leadFront + pcPbThickness;
          const G4double scintiBack = leadBack + pcScintiThickness;
          const auto suffix = std::to_string(layer);
          auto pcPbS = new G4Cons(
              "PCLeadS_" + suffix,
              leadFront * tanInner, leadFront * tanOuter,
              leadBack * tanInner, leadBack * tanOuter,
              pcPbThickness / 2., 0.*deg, 360.*deg);
          auto pcScintiS = new G4Cons(
              "PCScintiS_" + suffix,
              leadBack * tanInner, leadBack * tanOuter,
              scintiBack * tanInner, scintiBack * tanOuter,
              pcScintiThickness / 2., 0.*deg, 360.*deg);
          auto pcPbLV = new G4LogicalVolume(
              pcPbS, G4Material::GetMaterial("G4_Pb"),
              "PCLeadLV_" + suffix);
          auto pcScintiLV = new G4LogicalVolume(
              pcScintiS, ScintiMaterial, "PCScintiLV_" + suffix);
          pcPbLV->SetVisAttributes(pcPbVis);
          pcScintiLV->SetVisAttributes(pcScintiVis);
          new G4PVPlacement(
              0, G4ThreeVector(0., 0., 0.5 * (leadFront + leadBack)),
              pcPbLV, "PCLead", worldLV, false, layer, fCheckOverlaps);
          new G4PVPlacement(
              0, G4ThreeVector(0., 0., 0.5 * (leadBack + scintiBack)),
              pcScintiLV, "PCScinti", worldLV, false, layer,
              fCheckOverlaps);
        }
      }

      // Legacy PCUp means the -z endcap and is physically downstream (the
      // large BGOC opening).  Keep the names for ROOT compatibility.
      if (config.HasUpstreamPhotonCounter())
      {
        const G4double tanInner = std::tan(config.PCUpThetaInnerDeg() * deg);
        const G4double tanOuter = std::tan(config.PCUpThetaOuterDeg() * deg);
        for (G4int layer = 0; layer < pcNLayers; ++layer)
        {
          const G4double leadNear = pcZFront + layer * pcLayerThickness;
          const G4double leadFar = leadNear + pcPbThickness;
          const G4double scintiFar = leadFar + pcScintiThickness;
          const auto suffix = std::to_string(layer);
          // Local -z is the far face for a volume placed on global -z.
          auto pcUpPbS = new G4Cons(
              "PCUpLeadS_" + suffix,
              leadFar * tanInner, leadFar * tanOuter,
              leadNear * tanInner, leadNear * tanOuter,
              pcPbThickness / 2., 0.*deg, 360.*deg);
          auto pcUpScintiS = new G4Cons(
              "PCUpScintiS_" + suffix,
              scintiFar * tanInner, scintiFar * tanOuter,
              leadFar * tanInner, leadFar * tanOuter,
              pcScintiThickness / 2., 0.*deg, 360.*deg);
          auto pcUpPbLV = new G4LogicalVolume(
              pcUpPbS, G4Material::GetMaterial("G4_Pb"),
              "PCUpLeadLV_" + suffix);
          auto pcUpScintiLV = new G4LogicalVolume(
              pcUpScintiS, ScintiMaterial, "PCUpScintiLV_" + suffix);
          pcUpPbLV->SetVisAttributes(pcPbVis);
          pcUpScintiLV->SetVisAttributes(pcScintiVis);
          new G4PVPlacement(
              0, G4ThreeVector(0., 0., -0.5 * (leadNear + leadFar)),
              pcUpPbLV, "PCUpLead", worldLV, false, layer, fCheckOverlaps);
          new G4PVPlacement(
              0, G4ThreeVector(0., 0., -0.5 * (leadFar + scintiFar)),
              pcUpScintiLV, "PCUpScinti", worldLV, false, layer,
              fCheckOverlaps);
        }
      }
      G4cout << "PhotonCounterGeometrySummary mode=" << config.PhotonCounter()
             << " layers=" << pcNLayers
             << " z_cm=" << pcZFront / cm << "--" << pcZBack / cm
             << " down_theta_deg=" << config.PCDownThetaInnerDeg()
             << "--" << config.PCDownThetaOuterDeg()
             << " up_theta_deg=" << config.PCUpThetaInnerDeg()
             << "--" << config.PCUpThetaOuterDeg() << G4endl;
    }
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
  if (calorLV)
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

G4LogicalVolume *betaDetectorConstruction::ConstructCell(
    G4double rmin, G4double rmax,
    G4double thetaStart, G4double thetaSpan,
    G4double phiSpan, G4String mat)
{
  auto Material = G4Material::GetMaterial(mat);

  G4Sphere *cellS = new G4Sphere("Cell",
                                 rmin, rmax,
                                 0 * deg, phiSpan * deg,
                                 thetaStart * deg, thetaSpan * deg);

  cellLV = new G4LogicalVolume(
      cellS,    // its solid
      Material, // its material
      "Cell");  // its name

  return cellLV;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Transform3D betaDetectorConstruction::Rotate(G4double icrys, G4double phiSpan)
{
  G4ThreeVector pos = G4ThreeVector(0, 0, 0);
  G4RotationMatrix rotm = G4RotationMatrix();
  rotm.rotateZ(phiSpan * icrys * deg);
  G4Transform3D transform = G4Transform3D(rotm, pos);
  return transform;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void betaDetectorConstruction::ConstructBGOeggFrusta(
    G4String mat, std::vector<G4LogicalVolume *> &cells,
    G4LogicalVolume *mother, G4double zOffset)
{
  const auto &config = BetaConfig::Instance();
  const auto rings = BGOeggGeometry::BuildRings(config.NLayer());
  const G4double phiSpan = config.PhiWidthDeg();
  const G4double phiHalf = 0.5 * phiSpan * deg;
  auto material = G4Material::GetMaterial(mat);

  cells.assign(rings.size(), nullptr);
  G4int copyNo = 0;
  G4double activeVolume = 0.0;
  for (std::size_t ringIndex = 0; ringIndex < rings.size(); ++ringIndex)
  {
    const auto &ring = rings[ringIndex];
    const G4double thetaLow = ring.thetaLowDeg * deg;
    const G4double thetaHigh = ring.thetaHighDeg * deg;
    const G4double frontRadius = ring.frontRadiusMm * mm;
    const G4double rearRadius =
        (ring.frontRadiusMm + BGOeggGeometry::kCrystalLengthMm) * mm;

    auto vertex = [](G4double radius, G4double theta, G4double phi)
    {
      return G4ThreeVector(
          radius * std::sin(theta) * std::cos(phi),
          radius * std::sin(theta) * std::sin(phi),
          radius * std::cos(theta));
    };

    std::vector<G4ThreeVector> front = {
        vertex(frontRadius, thetaLow, -phiHalf),
        vertex(frontRadius, thetaLow, +phiHalf),
        vertex(frontRadius, thetaHigh, +phiHalf),
        vertex(frontRadius, thetaHigh, -phiHalf),
    };
    std::vector<G4ThreeVector> rear = {
        vertex(rearRadius, thetaLow, -phiHalf),
        vertex(rearRadius, thetaLow, +phiHalf),
        vertex(rearRadius, thetaHigh, +phiHalf),
        vertex(rearRadius, thetaHigh, -phiHalf),
    };

    G4ThreeVector normal =
        (front[1] - front[0]).cross(front[3] - front[0]).unit();
    const G4ThreeVector frontCenter =
        0.25 * (front[0] + front[1] + front[2] + front[3]);
    if (normal.dot(frontCenter) < 0.0)
      normal = -normal;

    const G4double frontPlane = normal.dot(front[0]);
    const G4double rearPlane = normal.dot(rear[0]);
    const G4double halfZ = 0.5 * (rearPlane - frontPlane);
    if (!(halfZ > 0.0))
      throw std::runtime_error("BGOegg frustum has non-positive thickness");

    const G4ThreeVector localX(0.0, 1.0, 0.0);
    const G4ThreeVector localY = normal.cross(localX).unit();
    const G4ThreeVector origin = normal * (0.5 * (frontPlane + rearPlane));

    std::vector<G4TwoVector> vertices;
    vertices.reserve(8);
    for (const auto &point : front)
    {
      const auto relative = point - origin;
      vertices.emplace_back(relative.dot(localX), relative.dot(localY));
    }
    for (const auto &point : rear)
    {
      const auto relative = point - origin;
      vertices.emplace_back(relative.dot(localX), relative.dot(localY));
    }

    const G4String solidName =
        "BGOeggCellSolid_" + std::to_string(ringIndex);
    auto solid = new G4GenericTrap(solidName, halfZ, vertices);
    activeVolume += solid->GetCubicVolume() * config.NSector();
    const G4String logicalName =
        "BGOeggCellLV_" + std::to_string(ringIndex);
    cells[ringIndex] = new G4LogicalVolume(solid, material, logicalName);

    const G4RotationMatrix baseRotation(localX, localY, normal);
    for (G4int sector = 0; sector < config.NSector(); ++sector)
    {
      G4RotationMatrix azimuth;
      // The local frustum spans -phi/2..+phi/2.  Offset by half a sector so
      // copy 0 keeps the historical global 0..phi sector convention.
      azimuth.rotateZ((sector + 0.5) * phiSpan * deg);
      const G4RotationMatrix rotation = azimuth * baseRotation;
      const G4ThreeVector position =
          azimuth * origin + G4ThreeVector(0.0, 0.0, zOffset);
      const G4Transform3D transform(rotation, position);
      fCellPV = new G4PVPlacement(
          transform, cells[ringIndex], "cellPhysical", mother, false,
          copyNo, fCheckOverlaps);
      ++copyNo;
    }
  }

  if (copyNo != config.NCells())
    throw std::runtime_error("BGOegg frustum copy count does not match configuration");
  G4cout << "BGOeggGeometrySummary model=" << config.GeometryModel()
         << " rings=" << rings.size()
         << " sectors=" << config.NSector()
         << " cells=" << copyNo
         << " theta_deg=" << rings.front().thetaLowDeg
         << "--" << rings.back().thetaHighDeg
         << " active_volume_liter=" << activeVolume / liter
         << " z_offset_cm=" << zOffset / cm << G4endl;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void betaDetectorConstruction::ConstructCalorimeter(
    G4double rmin, G4double rmax, G4String mat,
    std::vector<G4LogicalVolume *> &cells,
    G4LogicalVolume *mother)
{
  const auto &config = BetaConfig::Instance();
  const G4int nLayer = config.NLayer();
  const G4int nSector = config.NSector();
  const G4double phiSpan = config.PhiWidthDeg();
  const G4double activeThetaMin = config.ThetaMinDeg();
  const G4double activeThetaMax = config.ThetaMaxDeg();
  G4int CopyNo = 0;
  cells.assign(nLayer, nullptr);

  for (G4int j = 0; j < nLayer; j++)
  {
    G4double thetaStart = 0.0;
    G4double thetaStop = 0.0;
    if (config.EqualSolidAngle())
    {
      const G4double cosMin = std::cos(activeThetaMin * deg);
      const G4double cosMax = std::cos(activeThetaMax * deg);
      const G4double cosStart =
          cosMin + (cosMax - cosMin) * (static_cast<G4double>(j) / nLayer);
      const G4double cosStop =
          cosMin + (cosMax - cosMin) * (static_cast<G4double>(j + 1) / nLayer);
      thetaStart = std::acos(cosStart) / deg;
      thetaStop = std::acos(cosStop) / deg;
    }
    else
    {
      const G4double thetaSpan = (activeThetaMax - activeThetaMin) / nLayer;
      thetaStart = activeThetaMin + j * thetaSpan;
      thetaStop = thetaStart + thetaSpan;
    }
    cells[j] = ConstructCell(rmin, rmax, thetaStart,
                             thetaStop - thetaStart, phiSpan, mat);
    for (G4int icrys = 0; icrys < nSector; icrys++)
    {
      G4Transform3D transform = Rotate(icrys, phiSpan);
      fCellPV = new G4PVPlacement(transform, cells[j], "cellPhysical",
                                  mother, false, CopyNo, fCheckOverlaps);
      CopyNo++;
    }
  } // for j
}

void betaDetectorConstruction::ConstructSDandField()
{
  const G4int nLayer = BetaConfig::Instance().NLayer();

  // sensitive detectors -----------------------------------------------------
  auto sdManager = G4SDManager::GetSDMpointer();
  G4String SDname;

  auto calorimeter = new CalorimeterSD(SDname = "/calorimeter", "CalorimeterHC");
  sdManager->AddNewDetector(calorimeter);
  if (shapeFlag == 0)
  {
    cellLV->SetSensitiveDetector(calorimeter);
  }
  if (shapeFlag == 1)
  {
    auto target = new TargetSD(SDname = "/target", "TargetHC");
    sdManager->AddNewDetector(target);
    TargetLV->SetSensitiveDetector(target);

    auto th = new HodoscopeSD("/th", "THHC", false);
    sdManager->AddNewDetector(th);
    THsegLV->SetSensitiveDetector(th);

    auto tlc = new HodoscopeSD("/tlc", "TLCHC", true);
    sdManager->AddNewDetector(tlc);
    TLCsegLV->SetSensitiveDetector(tlc);

    for (G4int j = 0; j < nLayer; j++)
    {
      cellLV_Cal[j]->SetSensitiveDetector(calorimeter);
      // if (PhysicsFlag == 5)
      // {
      //   auto biasingOp =
      //       new betaBiasingOperator("CalorimeterBiasingOp", InelasticBias);
      //   biasingOp->AttachTo(cellLV_Cal[j]);
      // }
    }
  }

  // magnetic field ----------------------------------------------------------
  //  G4ThreeVector fieldValue;
  //  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  //  fMagFieldMessenger->SetVerboseLevel(1);
  //  Register the field messenger for deleting
  //  G4AutoDelete::Register(fMagFieldMessenger);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
