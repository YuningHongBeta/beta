#ifndef betaPrimaryGeneratorAction_h
#define betaPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4GenericMessenger;
class G4Event;

class betaPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  betaPrimaryGeneratorAction();
  ~betaPrimaryGeneratorAction() override;

  void GeneratePrimaries(G4Event* evt) override;

private:
  enum class Mode { kElectronLambda3Body, kPiMinus100MeV, kPiZero100MeV };

  void SetMode(const G4String& modeStr);
  G4String GetMode() const;

  void GenerateElectronLambda3Body(); // e-のエネルギー&方向をfGunにセット
  void GeneratePionFixedP(const G4String& pname); // π±/π0

private:
  G4ParticleGun* fGun = nullptr;

  // messenger
  G4GenericMessenger* fMessenger = nullptr;
  Mode fMode = Mode::kPiZero100MeV;

  // constants for e- spectrum (Lambda at rest, nu massless)
  G4double fMe = 0.0, fMp = 0.0, fMlam = 0.0, fMnu = 0.0;
  G4double fSmin = 0.0, fSmax = 0.0;
  G4double fWmax = 1.0;

  // fixed momentum for pions
  G4double fFixedP = 100.0; // will be multiplied by MeV in cc
};

#endif