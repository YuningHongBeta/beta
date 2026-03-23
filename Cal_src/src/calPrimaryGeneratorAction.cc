// calPrimaryGeneratorAction.cc
#include "calPrimaryGeneratorAction.hh"
#include "calConstant.hh"
// #include "ParamMan.hh"

#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4RandomDirection.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"

#include "TROOT.h"
#include "TH1.h"
#include "TTree.h"
#include "TFile.h"
#include "TRandom.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
calPrimaryGeneratorAction::calPrimaryGeneratorAction()
    : G4VUserPrimaryGeneratorAction(),
      fParticleGun(nullptr),
      Beamx(0.0), Beamy(0.0), Beamz(0.0),
      BeamPx(0.0), BeamPy(0.0), BeamPz(0.0),
      BeamMom(0.0),
      BeamE(0.0),
      decayFlag(0)
{
  G4int nofParticles = 1;
  fParticleGun = new G4ParticleGun(nofParticles);

  // default particle kinematic
  auto particleDefinition
      //= G4ParticleTable::GetParticleTable()->FindParticle("neutron");
      = G4ParticleTable::GetParticleTable()->FindParticle("e-");
  fParticleGun->SetParticleDefinition(particleDefinition);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
  fParticleGun->SetParticleEnergy(0. * MeV);   // init
  fParticleGun->SetParticleMomentum(0. * MeV); // init
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

calPrimaryGeneratorAction::~calPrimaryGeneratorAction()
{
  delete fParticleGun;
}

// allowed beta spectrum (no Fermi function): w(E) ∝ p E (E0 - E)^2
static G4double SampleElectronTotalEnergy_LambdaBeta(G4double Mlam,
                                                     G4double Mp,
                                                     G4double me)
{
  // Endpoint for electron total energy when neutrino is soft:
  // E_e,max = (Mlam^2 + me^2 - Mp^2) / (2 Mlam)
  const G4double Ee_min = me;
  const G4double Ee_max = (Mlam * Mlam + me * me - Mp * Mp) / (2.0 * Mlam);

  if (Ee_max <= Ee_min)
    return Ee_min; // kinematically impossible case

  // Precompute wmax by scan (robust, avoids arbitrary normalization constants)
  G4double wmax = 0.0;
  const int Nscan = 20000;
  for (int i = 0; i <= Nscan; ++i)
  {
    const G4double Ee = Ee_min + (Ee_max - Ee_min) * (G4double(i) / Nscan);
    const G4double pe2 = std::max(0.0, Ee * Ee - me * me);
    const G4double pe = std::sqrt(pe2);
    const G4double w = pe * Ee * std::pow(Ee_max - Ee, 2);
    wmax = std::max(wmax, w);
  }
  wmax *= 1.2; // safety margin

  // Rejection sampling
  while (true)
  {
    const G4double Ee = Ee_min + (Ee_max - Ee_min) * G4UniformRand();
    const G4double pe2 = std::max(0.0, Ee * Ee - me * me);
    const G4double pe = std::sqrt(pe2);
    const G4double w = pe * Ee * std::pow(Ee_max - Ee, 2);

    if (G4UniformRand() * wmax < w)
      return Ee;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void calPrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent)
{
  auto particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("e-");

  ////////////////////////////////////////
  // Decay flag
  ////////////////////////////////////////
  G4double rnd = gRandom->Uniform(0, 1);
  if (rnd < 0.5)
  {
    decayFlag = 1;
  }
  else if (rnd < 1.0)
  {
    decayFlag = 2;
  }
  decayFlag = 0;

  if (decayFlag == 0)
  {
    ////////////////////////////////////////
    // 100 MeV electron
    ////////////////////////////////////////
    particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("pi-");
    // particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("e-");
    // BeamE = 100*MeV;
    BeamMom = 100. * MeV;
  }

  if (decayFlag == 1)
  {
    ////////////////////////////////////////
    // Lambda strong decay
    ////////////////////////////////////////
    particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("lambda");
    BeamE = 0 * MeV;
    BeamMom = 100. * MeV;
  }
  else if (decayFlag == 2)
  {
    ////////////////////////////////////////
    // Lambda beta decay
    ////////////////////////////////////////
    particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("e-");
    const G4double Ee = SampleElectronTotalEnergy_LambdaBeta(M_Lambda, M_proton, M_e);
    const G4double Eke = std::max(0.0, Ee - M_e);
    BeamE = Eke * MeV;
    BeamMom = std::sqrt(Ee * Ee - M_e * M_e) * MeV;
  }

  if (shapeFlag == 0)
  {                                                         // for Box shape
    beamPos = G4ThreeVector(0., 0., -(absoThickness / 2.)); // from Box edge
    direction = G4ThreeVector(0., 0., 1.);                  // z direction
  }
  if (shapeFlag == 1)
  {                                      // for Shere shape
    beamPos = G4ThreeVector(0., 0., 0.); // world center
    direction = G4RandomDirection();     // uniform from (0,0,0)
  }

  // direction = G4ThreeVector(0.,1.,1.); // z direction

  fParticleGun->SetParticleDefinition(particleDefinition);
  // fParticleGun->SetParticleEnergy(BeamE);
  fParticleGun->SetParticlePosition(beamPos);
  fParticleGun->SetParticleMomentum(BeamMom);
  fParticleGun->SetParticleMomentumDirection(direction);
  fParticleGun->GeneratePrimaryVertex(anEvent);

  ////////////////////////
  // save primary information
  ////////////////////////
  TTree *tree = (TTree *)gROOT->FindObject("tree");
  tree->SetBranchAddress("beamE", &BeamE);
  tree->SetBranchAddress("decayFlag", &decayFlag);

  G4double worldZHalfLength = 0.;
  auto worldLV = G4LogicalVolumeStore::GetInstance()->GetVolume("World");

  // Check that the world volume has box shape
  G4Box *worldBox = nullptr;
  if (worldLV)
  {
    worldBox = dynamic_cast<G4Box *>(worldLV->GetSolid());
  }

  if (worldBox)
  {
    worldZHalfLength = worldBox->GetZHalfLength();
  }
  else
  {
    G4ExceptionDescription msg;
    msg << "World volume of box shape not found." << G4endl;
    msg << "Perhaps you have changed geometry." << G4endl;
    msg << "The gun will be place in the center.";
    G4Exception("calPrimaryGeneratorAction::GeneratePrimaries()",
                "MyCode0002", JustWarning, msg);
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void calPrimaryGeneratorAction::init()
{
  m12 = 0;
  m23 = 0;
  Ee = 0.;
}
G4double calPrimaryGeneratorAction::incident(double ma, double pa, double ua, double va,
                                             double mb, double pb, double ub, double vb)
{
  // incident
  pa_z = pa / sqrt(1 + ua * ua + va * va);
  pa_x = pa_z * ua;
  pa_y = pa_z * va;
  pA.setX(pa_x);
  pA.setY(pa_y);
  pA.setZ(pa_z);
  Ea = sqrt(ma * ma + pA.mag() * pA.mag());

  pb_z = pb / sqrt(1 + ub * ub + vb * vb);
  pb_x = pb_z * ub;
  pb_y = pb_z * vb;
  pB.setX(pb_x);
  pB.setY(pb_y);
  pB.setZ(pb_z);
  Eb = sqrt(mb * mb + pB.mag() * pB.mag());

  // CM system
  psum = pA + pB;
  psum1 = psum.mag();
  Esum1 = (Ea + Eb);
  beta.setX(psum.getX() / Esum1);
  beta.setY(psum.getY() / Esum1);
  beta.setZ(psum.getZ() / Esum1);
  // gamma = 1/sqrt(1-beta.Mag()*beta.Mag());
  Ecm = sqrt(Esum1 * Esum1 - psum1 * psum1);
  return Ecm;
}

G4double calPrimaryGeneratorAction::ThreeBodyDecay(G4double m1, G4double m2, G4double m3, G4double Ecm)
{
  M = Ecm;
  m1m2 = m1 + m2; // m1 + m2
  Mm3 = Ecm - m3; // M - m3

  m12 = gRandom->Uniform(m1m2, Mm3);
  E2 = (m12 * m12 - m1 * m1 + m2 * m2) / (2 * m12);
  E3 = (M * M - m12 * m12 - m3 * m3) / (2 * m12);
  tmp1 = sqrt(E2 * E2 - m2 * m2) - sqrt(E3 * E3 - m3 * m3);
  tmp2 = sqrt(E2 * E2 - m2 * m2) + sqrt(E3 * E3 - m3 * m3);
  m23max = pow((E2 + E3), 2) - pow(tmp1, 2);
  m23min = pow((E2 + E3), 2) - pow(tmp2, 2);

  m23 = gRandom->Uniform(sqrt(m23min), sqrt(m23max));
  tmp3 = (M * M - (m12 + m3) * (m12 + m3)) * (M * M - (m12 - m3) * (m12 - m3));
  p3 = sqrt(tmp3) / (2 * M);      // p3
  E3cm = sqrt(m3 * m3 + p3 * p3); // E3cm

  return E3cm;
}
