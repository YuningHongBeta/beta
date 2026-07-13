#ifndef HodoscopeSD_h
#define HodoscopeSD_h 1

#include "G4VSensitiveDetector.hh"
#include "HodoscopeHit.hh"
#include "globals.hh"

#include <vector>

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

class HodoscopeSD : public G4VSensitiveDetector
{
public:
  HodoscopeSD(const G4String &name, const G4String &hitsCollectionName,
              G4bool calculateCherenkov);
  ~HodoscopeSD() override = default;

  void Initialize(G4HCofThisEvent *hce) override;
  G4bool ProcessHits(G4Step *step, G4TouchableHistory *history) override;

private:
  HodoscopeHitsCollection *fHitsCollection = nullptr;
  G4int fHCID = -1;
  G4bool fCalculateCherenkov = false;
  std::vector<HodoscopeHit *> fHitByCopyNo;
};

#endif