// CalorimeterHit.hh

#ifndef CalorimeterHit_h
#define CalorimeterHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

class G4AttDef;
class G4AttValue;

class CalorimeterHit : public G4VHit
{
  public:
    CalorimeterHit();
    CalorimeterHit(G4int cellID);
    CalorimeterHit(const CalorimeterHit &right);
    virtual ~CalorimeterHit();

    const CalorimeterHit& operator=(const CalorimeterHit &right);
    G4bool operator==(const CalorimeterHit &right) const;
    
    inline void *operator new(size_t);
    inline void operator delete(void *aHit);
    
    virtual void Draw();
    virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
    virtual std::vector<G4AttValue>* CreateAttValues() const;
    virtual void Print();
    
    void SetCellID(G4int z) { fCellID = z; }
    G4int GetCellID() const { return fCellID; }

    void SetEdep(G4double de) { fEdep = de; }
    void AddEdep(G4double de) { fEdep += de; }
    G4double GetEdep() const { return fEdep; }

    void SetPos(G4ThreeVector xyz) { fPos = xyz; }
    G4ThreeVector GetPos() const { return fPos; }

    void SetRot(G4RotationMatrix rmat) { fRot = rmat; }
    G4RotationMatrix GetRot() const { return fRot; }

    void SetLogV(G4LogicalVolume* val) { fPLogV = val; }
    const G4LogicalVolume* GetLogV() const { return fPLogV; }
    
  private:
    G4int fCellID;
    G4double fEdep;
    G4ThreeVector fPos;
    G4RotationMatrix fRot;
    const G4LogicalVolume* fPLogV;
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

inline void CalorimeterHit::operator delete(void* aHit)
{
  CalorimeterHitAllocator->FreeSingle((CalorimeterHit*) aHit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
