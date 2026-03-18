#ifndef betaEventAction_hh
#define betaEventAction_hh 1

#include "G4UserEventAction.hh"
#include "globals.hh"

#include <memory>
#include <vector>

class G4Event;
class G4GenericMessenger;

class betaEventAction : public G4UserEventAction
{
public:
  betaEventAction();
  ~betaEventAction() override;

  void BeginOfEventAction(const G4Event *evt) override;
  void EndOfEventAction(const G4Event*) override;

  // ---- SteppingAction からの加算用（コンパイルエラー回避 + 使うならここに足す）----
  void AddEdepCell(double edep)   { fEvtEdepCell   += edep; } // edep: internal unit
  void AddEdepTarget(double edep) { fEvtEdepTarget += edep; } // edep: internal unit
  void SetWatchEvtID(int evtID)   { fWatchEvtID = evtID; }

  // ---- calarr 用（RunAction が vector column を作る時に参照）----
  const std::vector<double>& GetDE_MeV() const { return *fDE_MeV; }
  const std::vector<int>&    GetPID()   const { return *fPID;   }

private:
  double fEvtEdepCell   = 0.0; // internal unit
  double fEvtEdepTarget = 0.0; // internal unit
  int fWatchEvtID = -1;

  std::vector<double>* fDE_MeV; // size=225, MeV, heap allocated
  std::vector<int>*    fPID;    // size=225, PDG code, heap allocated

  std::unique_ptr<G4GenericMessenger> fMessenger;

  G4bool   fDoTargetTrace     = false; // ターゲット内の寄与内訳を毎回printするか
  G4bool   fSaveRndmOnPeak    = false; // ピーク窓に入ったら rndm を保存するか
  G4double fPeakEmin_MeV      = 6.0;   // MeV
  G4double fPeakEmax_MeV      = 10.0;   // MeV
  G4String fRndmDir           = ".";   // 保存先
  G4String fRndmPrefix        = "rndm_evt"; // ファイル名prefix
};

#endif