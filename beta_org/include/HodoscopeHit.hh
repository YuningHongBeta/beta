#ifndef HodoscopeHit_h
#define HodoscopeHit_h 1

#include "G4THitsCollection.hh"
#include "G4VHit.hh"
#include "globals.hh"

class HodoscopeHit : public G4VHit
{
public:
  HodoscopeHit() = default;
  ~HodoscopeHit() override = default;

  void SetCopyNo(G4int value) { fCopyNo = value; }
  void SetTime(G4double value) { fTime = value; }
  void SetPDG(G4int value) { fPDG = value; }
  void UpdateCherenkovTime(G4double value)
  {
    if (value < fCherenkovTime)
      fCherenkovTime = value;
  }
  void AddEdep(G4double value) { fEdep += value; }
  void AddChargedPath(G4double value) { fChargedPath += value; }
  void AddCherenkovPath(G4double value) { fCherenkovPath += value; }
  void AddCherenkovExpectedPhotons(G4double value) { fCherenkovExpectedPhotons += value; }

  G4int GetCopyNo() const { return fCopyNo; }
  G4double GetTime() const { return fTime; }
  G4int GetPDG() const { return fPDG; }
  G4double GetEdep() const { return fEdep; }
  G4double GetChargedPath() const { return fChargedPath; }
  G4double GetCherenkovPath() const { return fCherenkovPath; }
  G4double GetCherenkovTime() const { return fCherenkovTime; }
  G4double GetCherenkovExpectedPhotons() const { return fCherenkovExpectedPhotons; }

private:
  G4int fCopyNo = -1;
  G4double fTime = 0.0;
  G4int fPDG = 0;
  G4double fEdep = 0.0;
  G4double fChargedPath = 0.0;
  G4double fCherenkovPath = 0.0;
  G4double fCherenkovTime = 1.0e30;
  G4double fCherenkovExpectedPhotons = 0.0;
};

using HodoscopeHitsCollection = G4THitsCollection<HodoscopeHit>;

#endif