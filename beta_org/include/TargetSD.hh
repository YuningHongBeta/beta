#ifndef TargetSD_h
#define TargetSD_h 1

#include "G4VSensitiveDetector.hh"
#include "globals.hh"

// TargetHit.hh の中で TargetHit と TargetHitsCollection が定義されている前提
#include "TargetHit.hh"

#include <map>
#include <string>

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

class TargetSD : public G4VSensitiveDetector
{
public:
  TargetSD(const G4String& name, const G4String& hitsCollectionName);
  ~TargetSD() override = default;

  void Initialize(G4HCofThisEvent* hce) override;
  G4bool ProcessHits(G4Step* step, G4TouchableHistory* history) override;
  void EndOfEvent(G4HCofThisEvent* hce) override;

  // 9MeVピーク窓（MeV）。必要ならmacro的に変えられるように。
  void SetPeakWindowMeV(double lo, double hi) { fPeakLoMeV = lo; fPeakHiMeV = hi; }

private:
  TargetHitsCollection* fHitsCollection = nullptr;
  G4int fHCID = -1;

  // 1 event = 1 hit
  TargetHit* fHit = nullptr;

  // --- trace accumulators (internal unit: energy) ---
  G4double fEdepTot = 0.0;
  std::map<G4int, G4double> fEdepByPDG;
  std::map<std::string, G4double> fEdepByProc;

  // peak window in MeV
  G4double fPeakLoMeV = 8.5;
  G4double fPeakHiMeV = 9.5;
};

#endif