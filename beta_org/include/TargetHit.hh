#ifndef TargetHit_h
#define TargetHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

class TargetHit : public G4VHit
{
public:
  TargetHit() = default;
  ~TargetHit() override = default;

  // --- setters ---
  void SetTime(G4double v) { fTime = v; }
  void SetEdep(G4double v) { fEdep = v; }
  void SetPDG(G4int v)     { fPDG  = v; }
  void SetMomentum(const G4ThreeVector& v) { fMom = v; }

  // --- getters ---
  G4double GetTime() const { return fTime; }
  G4double GetEdep() const { return fEdep; }
  G4int    GetPDG()  const { return fPDG; }
  const G4ThreeVector& GetMomentum() const { return fMom; }

  // (optional) add edep
  void AddEdep(G4double v) { fEdep += v; }

private:
  G4double      fTime = 0.0;
  G4double      fEdep = 0.0;
  G4int         fPDG  = 0;
  G4ThreeVector fMom  = G4ThreeVector(0,0,0);
};

// hits collection typedef
using TargetHitsCollection = G4THitsCollection<TargetHit>;

#endif