#include "betaPrimaryGeneratorAction.hh"
#include "BetaConfig.hh"
#include "Constant.hh"

#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
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

  // ---- set up messenger: clean, beam-only, or signal+beam modes
  fMessenger = new G4GenericMessenger(this, "/gen/", "Primary generator control");
  fMessenger->DeclareMethod("mode", &betaPrimaryGeneratorAction::SetMode)
    .SetGuidance("Select primary: e|pim|pi0|beam|e_beam|pim_beam|pi0_beam")
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

  SetMode(BetaConfig::Instance().Primary());
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
  fOverlayBeam = false;

  if (m == "e" || m == "electron") {
    fMode = Mode::kElectronLambda3Body;
  } else if (m == "pim" || m == "pi-" || m == "piminus") {
    fMode = Mode::kPiMinus100MeV;
  } else if (m == "pi0") {
    fMode = Mode::kPiZero100MeV;
  } else if (m == "beam") {
    fMode = Mode::kBeamOnly;
    fOverlayBeam = true;
  } else if (m == "e_beam") {
    fMode = Mode::kElectronLambda3Body;
    fOverlayBeam = true;
  } else if (m == "pim_beam") {
    fMode = Mode::kPiMinus100MeV;
    fOverlayBeam = true;
  } else if (m == "pi0_beam") {
    fMode = Mode::kPiZero100MeV;
    fOverlayBeam = true;
  } else {
    G4Exception("betaPrimaryGeneratorAction::SetMode", "Gen001", JustWarning,
                ("Unknown mode: " + modeStr +
                 " (use e|pim|pi0|beam|e_beam|pim_beam|pi0_beam)").c_str());
    fMode = Mode::kElectronLambda3Body;
  }
}

G4String betaPrimaryGeneratorAction::GetMode() const
{
  switch (fMode) {
    case Mode::kElectronLambda3Body: return "e";
    case Mode::kPiMinus100MeV:       return "pim";
    case Mode::kPiZero100MeV:        return "pi0";
    case Mode::kBeamOnly:            return "beam";
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

  const auto &config = BetaConfig::Instance();
  const G4double p = (pname == "pi-" ? config.PiMinusMomentumMeVC()
                                         : config.PiZeroMomentumMeVC()) * MeV;
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

void betaPrimaryGeneratorAction::GenerateIncidentBeam(G4Event* evt)
{
  const auto &config = BetaConfig::Instance();
  auto* part = G4ParticleTable::GetParticleTable()->FindParticle(
      config.BeamParticle());
  if (!part) {
    G4Exception("betaPrimaryGeneratorAction::GenerateIncidentBeam", "Gen003",
                FatalException,
                ("Beam particle not found: " + config.BeamParticle()).c_str());
    return;
  }

  const G4double momentum = config.BeamMomentumMeVC() * MeV;
  const G4double mass = part->GetPDGMass();
  const G4double totalEnergy = std::sqrt(momentum * momentum + mass * mass);
  const G4double beta = momentum / totalEnergy;
  const G4double upstreamDistance =
      config.TargetLengthMm() * mm / 2 + 0.1 * mm;

  fGun->SetParticleDefinition(part);
  fGun->SetParticleEnergy(totalEnergy - mass);
  // In this geometry theta_min is the upstream/beam-entrance opening (+z).
  fGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., -1.));
  fGun->SetParticlePosition(G4ThreeVector(0., 0., upstreamDistance));
  // Define t=0 at the target center, where the clean signal vertex is placed.
  fGun->SetParticleTime(-upstreamDistance / (beta * c_light));
  fGun->GeneratePrimaryVertex(evt);
}

void betaPrimaryGeneratorAction::GeneratePrimaries(G4Event* evt)
{
  if (fMode != Mode::kBeamOnly) {
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
      case Mode::kBeamOnly:
        break;
    }

    fGun->SetParticlePosition(G4ThreeVector(0,0,0));
    fGun->SetParticleTime(0.0);
    fGun->GeneratePrimaryVertex(evt);
  }

  if (fOverlayBeam)
    GenerateIncidentBeam(evt);
}
