#ifndef CalorimeterSD_hh
#define CalorimeterSD_hh 1

#include "G4VSensitiveDetector.hh"
#include "CalorimeterHit.hh"
#include <unordered_set>
#include <cstdint>

#include <vector>
#include <unordered_map>

class G4Step;
class G4HCofThisEvent;

class CalorimeterSD : public G4VSensitiveDetector
{
public:
  CalorimeterSD(const G4String& name, const G4String& hitsCollectionName);
  ~CalorimeterSD() override = default;

  void Initialize(G4HCofThisEvent* hce) override;
  G4bool ProcessHits(G4Step* step, G4TouchableHistory*) override;

  // void SetTraceEnter(bool v) { fTraceEnter = v; }

private:
  CalorimeterHitsCollection* fHitsCollection = nullptr;
  G4int fHCID = -1;

  // 1 segment(copyNo) -> 1 hit
  std::vector<CalorimeterHit*> fHitByCopyNo;

  // ★ copyNo ごとに「pdg -> edep(内部単位)」を積算して最大寄与pdgを決める
  std::vector<std::unordered_map<int, G4double>> fEdepByPdgByCopyNo;
  std::vector<G4double> fMaxEdepByCopyNo;
  std::vector<int>      fMaxPdgByCopyNo;

  void EnsureSize(int copyNo);

  // bool fTraceEnter = false; // 侵入ログを出すかどうか

  // std::unordered_set<std::uint64_t> fPrinted;
};

#endif