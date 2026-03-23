#pragma once

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4String.hh"

class CalorimeterHit : public G4VHit {
public:
  CalorimeterHit() = default;
  ~CalorimeterHit() override = default;

  CalorimeterHit(const CalorimeterHit&) = default;
  CalorimeterHit& operator=(const CalorimeterHit&) = default;

  inline void* operator new(size_t);
  inline void  operator delete(void*);

  // Setters
  void SetCopyNo(G4int v)             { fCopyNo = v; }
  void SetTime(G4double v)            { fTime = v; }          // global time
  void SetEdep(G4double v)            { fEdep = v; }
  void SetPDG(G4int v)                { fPDG = v; }
  void SetMomentum(const G4ThreeVector& v) { fP = v; }
  void SetCreator(const G4String& s)  { fCreator = s; }
  void SetOriginType(G4int v)         { fOriginType = v; }

  void AddEdep(G4double v) { fEdep += v; }

  // Getters
  G4int GetCopyNo() const             { return fCopyNo; }
  G4double GetTime() const            { return fTime; }
  G4double GetEdep() const            { return fEdep; }
  G4int GetPDG() const                { return fPDG; }
  const G4ThreeVector& GetMomentum() const { return fP; }
  const G4String& GetCreator() const  { return fCreator; }
  G4int GetOriginType() const         { return fOriginType; }

private:
  G4int        fCopyNo = -1;
  G4double     fTime   = 0.0;
  G4double     fEdep   = 0.0;
  G4int        fPDG    = 0;
  G4ThreeVector fP     = G4ThreeVector();
  G4String     fCreator = "Unknown";
  G4int        fOriginType = 0; // 0:unknown, 1:evap, 2:fission
};

using CalorimeterHitsCollection = G4THitsCollection<CalorimeterHit>;

extern G4ThreadLocal G4Allocator<CalorimeterHit>* CalorimeterHitAllocator;

inline void* CalorimeterHit::operator new(size_t)
{
  if (!CalorimeterHitAllocator) {
    CalorimeterHitAllocator = new G4Allocator<CalorimeterHit>;
  }
  return (void*)CalorimeterHitAllocator->MallocSingle();
}

inline void CalorimeterHit::operator delete(void* hit)
{
  CalorimeterHitAllocator->FreeSingle((CalorimeterHit*)hit);
}