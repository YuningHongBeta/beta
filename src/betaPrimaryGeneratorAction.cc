#include "betaPrimaryGeneratorAction.hh"

#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4Event.hh"
#include "G4ThreeVector.hh"
#include "G4RandomDirection.hh"
#include "G4GenericMessenger.hh"
#include "Randomize.hh"
#include <cmath>
#include <algorithm>

namespace {

// Källén function: λ(a,b,c) = a^2 + b^2 + c^2 - 2ab - 2ac - 2bc
inline G4double Kallen(G4double a, G4double b, G4double c)
{
  return a*a + b*b + c*c - 2.0*(a*b + a*c + b*c);
}

// Lorentz boost of (E, p) by velocity beta (3-vector)
inline void Boost(const G4ThreeVector& beta,
                  G4double Ein, const G4ThreeVector& pin,
                  G4double& Eout, G4ThreeVector& pout)
{
  const G4double b2 = beta.mag2();
  if (b2 < 1e-18) { Eout = Ein; pout = pin; return; }
  const G4double gamma = 1.0 / std::sqrt(1.0 - b2);
  const G4double bp    = beta.dot(pin);

  Eout = gamma * (Ein + bp);
  pout = pin + (((gamma - 1.0) * bp / b2) + gamma * Ein) * beta;
}

} // namespace

betaPrimaryGeneratorAction::betaPrimaryGeneratorAction()
{
  fGun = new G4ParticleGun(1);
  fGun->SetParticlePosition(G4ThreeVector(0,0,0));
  fGun->SetParticleMomentumDirection(G4RandomDirection());

  // ---- set up messenger: /gen/mode [e|pim|pi0]
  fMessenger = new G4GenericMessenger(this, "/gen/", "Primary generator control");
  fMessenger->DeclareMethod("mode", &betaPrimaryGeneratorAction::SetMode)
    .SetGuidance("Select primary: e (Lambda->p e nu 3body), pim (pi- 100MeV/c), pi0 (pi0 100MeV/c)")
    .SetParameterName("mode", false);

  // ---- masses for e-spectrum
  auto* table = G4ParticleTable::GetParticleTable();

  auto* ePart = table->FindParticle("e-");
  fMe = ePart->GetPDGMass();

  auto* pPart = table->FindParticle("proton");
  fMp = pPart ? pPart->GetPDGMass() : 938.2720813 * MeV;

  auto* lamPart = table->FindParticle("lambda");
  fMlam = lamPart ? lamPart->GetPDGMass() : 1115.683 * MeV;

  fMnu = 0.0;

  fSmin = (fMe + fMnu) * (fMe + fMnu);
  fSmax = (fMlam - fMp) * (fMlam - fMp);

  // Precompute Wmax for rejection sampling (simple scan)
  fWmax = 0.0;
  const int Nscan = 20000;
  const G4double M2  = fMlam * fMlam;
  const G4double mp2 = fMp * fMp;
  const G4double me2 = fMe * fMe;

  for (int i = 0; i <= Nscan; ++i) {
    const G4double s = fSmin + (fSmax - fSmin) * (G4double(i) / Nscan);
    const G4double l1 = Kallen(M2, mp2, s);
    const G4double l2 = Kallen(s, me2, 0.0);
    if (l1 <= 0.0 || l2 <= 0.0 || s <= 0.0) continue;
    const G4double w = std::sqrt(l1) * std::sqrt(l2) / std::sqrt(s);
    if (w > fWmax) fWmax = w;
  }
  fWmax *= 1.2; // safety margin

  // default
  SetMode("e");
}

betaPrimaryGeneratorAction::~betaPrimaryGeneratorAction()
{
  delete fMessenger;
  delete fGun;
}

void betaPrimaryGeneratorAction::SetMode(const G4String& modeStr)
{
  G4String m = modeStr;
  m.toLower();

  if (m == "e" || m == "electron") {
    fMode = Mode::kElectronLambda3Body;
  } else if (m == "pim" || m == "pi-" || m == "piminus") {
    fMode = Mode::kPiMinus100MeV;
  } else if (m == "pi0") {
    fMode = Mode::kPiZero100MeV;
  } else {
    G4Exception("betaPrimaryGeneratorAction::SetMode", "Gen001", JustWarning,
                ("Unknown mode: " + modeStr + " (use e|pim|pi0)").c_str());
    fMode = Mode::kElectronLambda3Body;
  }
}

G4String betaPrimaryGeneratorAction::GetMode() const
{
  switch (fMode) {
    case Mode::kElectronLambda3Body: return "e";
    case Mode::kPiMinus100MeV:       return "pim";
    case Mode::kPiZero100MeV:        return "pi0";
  }
  return "e";
}

void betaPrimaryGeneratorAction::GeneratePionFixedP(const G4String& pname)
{
  auto* table = G4ParticleTable::GetParticleTable();
  auto* part  = table->FindParticle(pname);
  if (!part) {
    G4Exception("betaPrimaryGeneratorAction::GeneratePionFixedP", "Gen002", FatalException,
                ("Particle not found: " + pname).c_str());
    return;
  }

  fGun->SetParticleDefinition(part);

  const G4double p = fFixedP * MeV;          // 100 MeV/c を c=1で扱う
  const G4double m = part->GetPDGMass();
  const G4double E = std::sqrt(p*p + m*m);  // total energy
  const G4double Ek = E - m;

  fGun->SetParticleEnergy(std::max(0.0, Ek));
  fGun->SetParticleMomentumDirection(G4RandomDirection());
}

void betaPrimaryGeneratorAction::GenerateElectronLambda3Body()
{
  auto* table = G4ParticleTable::GetParticleTable();
  auto* ePart = table->FindParticle("e-");
  fGun->SetParticleDefinition(ePart);

  const G4double M2  = fMlam * fMlam;
  const G4double mp2 = fMp * fMp;
  const G4double me2 = fMe * fMe;

  // sample s = m^2_{eν}
  G4double s = 0.0;
  while (true) {
    s = fSmin + (fSmax - fSmin) * G4UniformRand();
    const G4double l1 = Kallen(M2, mp2, s);
    const G4double l2 = Kallen(s, me2, 0.0);
    if (l1 <= 0.0 || l2 <= 0.0 || s <= 0.0) continue;
    const G4double w = std::sqrt(l1) * std::sqrt(l2) / std::sqrt(s);
    if (G4UniformRand() * fWmax < w) break;
  }

  // proton direction in Lambda rest
  const G4ThreeVector n_p = G4RandomDirection();

  // proton momentum magnitude in Lambda rest
  const G4double p_p = std::sqrt(std::max(0.0, Kallen(M2, mp2, s))) / (2.0 * fMlam);
  const G4double E_p = std::sqrt(mp2 + p_p*p_p);

  // eν system in Lambda rest
  const G4double E_sys = fMlam - E_p;
  const G4ThreeVector p_sys = -p_p * n_p;
  const G4ThreeVector beta_sys = p_sys / E_sys;

  // electron in eν rest
  const G4double sqrt_s = std::sqrt(s);
  const G4double p_e_star =
      std::sqrt(std::max(0.0, Kallen(s, me2, 0.0))) / (2.0 * sqrt_s);
  const G4double E_e_star = std::sqrt(me2 + p_e_star*p_e_star);

  const G4ThreeVector n_star = G4RandomDirection();
  const G4ThreeVector p_e_star_vec = p_e_star * n_star;

  // boost to Lambda rest
  G4double E_e = 0.0;
  G4ThreeVector p_e;
  Boost(beta_sys, E_e_star, p_e_star_vec, E_e, p_e);

  const G4double Ek = E_e - fMe;
  fGun->SetParticleEnergy(std::max(0.0, Ek));
  fGun->SetParticleMomentumDirection(p_e.unit());
}

void betaPrimaryGeneratorAction::GeneratePrimaries(G4Event* evt)
{
  switch (fMode) {
    case Mode::kElectronLambda3Body:
      GenerateElectronLambda3Body();
      break;
    case Mode::kPiMinus100MeV:
      GeneratePionFixedP("pi-");
      break;
    case Mode::kPiZero100MeV:
      GeneratePionFixedP("pi0");
      break;
  }

  fGun->SetParticlePosition(G4ThreeVector(0,0,0));
  fGun->GeneratePrimaryVertex(evt);
}
