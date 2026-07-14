#ifndef betaEventAction_hh
#define betaEventAction_hh 1

#include "G4UserEventAction.hh"
#include "globals.hh"

#include <memory>
#include <unordered_set>
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
  void AddEdepPCDown(double edep)
  {
    fEvtEdepPC += edep;
    fEvtEdepPCDown += edep;
  }
  void AddEdepPCUp(double edep)
  {
    fEvtEdepPC += edep;
    fEvtEdepPCUp += edep;
  }
  void RecordPCGammaEntrance(bool upstream, int trackID, double energy);
  void SetWatchEvtID(int evtID)   { fWatchEvtID = evtID; }

  // ---- calarr 用（RunAction が vector column を作る時に参照）----
  std::vector<double>& GetDE_MeV() { return fDE_MeV; }
  std::vector<int>&    GetPID()    { return fPID; }
  std::vector<double>& GetTHDE_MeV() { return fTHDE_MeV; }
  std::vector<double>& GetTHTime_ns() { return fTHTime_ns; }
  std::vector<double>& GetTHTimeLeft_ns() { return fTHTimeLeft_ns; }
  std::vector<double>& GetTHTimeRight_ns() { return fTHTimeRight_ns; }
  std::vector<double>& GetTHTimeLeftMinusRight_ns() { return fTHTimeLeftMinusRight_ns; }
  std::vector<double>& GetTHZReco_mm() { return fTHZReco_mm; }
  std::vector<double>& GetTHPathTruth_mm() { return fTHPathTruth_mm; }
  std::vector<double>& GetTLCDE_MeV() { return fTLCDE_MeV; }
  std::vector<double>& GetTLCCherenkovTime_ns() { return fTLCCherenkovTime_ns; }
  std::vector<double>& GetTLCPath_mm() { return fTLCPath_mm; }
  std::vector<double>& GetTLCCherenkovPath_mm() { return fTLCCherenkovPath_mm; }
  std::vector<double>& GetTLCCherenkovExpectedPhotons() { return fTLCCherenkovExpectedPhotons; }

private:
  double fEvtEdepCell   = 0.0; // internal unit
  double fEvtEdepTarget = 0.0; // internal unit
  double fEvtEdepPC     = 0.0; // enabled photon-counter plastic layers
  double fEvtEdepPCDown = 0.0;
  double fEvtEdepPCUp   = 0.0;
  int fPCGammaN = 0;
  int fPCGammaDownN = 0;
  int fPCGammaUpN = 0;
  double fPCGammaEnergy = 0.0;
  double fPCGammaDownEnergy = 0.0;
  double fPCGammaUpEnergy = 0.0;
  double fPCGammaMaxEnergy = 0.0;
  std::unordered_set<int> fPCGammaDownTrackIDs;
  std::unordered_set<int> fPCGammaUpTrackIDs;
  int fWatchEvtID = -1;

  std::vector<double> fDE_MeV;
  std::vector<int>    fPID;
  std::vector<double> fTHDE_MeV;
  std::vector<double> fTHTime_ns;
  std::vector<double> fTHTimeLeft_ns;
  std::vector<double> fTHTimeRight_ns;
  std::vector<double> fTHTimeLeftMinusRight_ns;
  std::vector<double> fTHZReco_mm;
  std::vector<double> fTHPathTruth_mm;
  std::vector<double> fTLCDE_MeV;
  std::vector<double> fTLCCherenkovTime_ns;
  std::vector<double> fTLCPath_mm;
  std::vector<double> fTLCCherenkovPath_mm;
  std::vector<double> fTLCCherenkovExpectedPhotons;

  G4int fCalHCID = -1;
  G4int fTargetHCID = -1;
  G4int fTHHCID = -1;
  G4int fTLCHCID = -1;

  std::unique_ptr<G4GenericMessenger> fMessenger;

  G4bool   fDoTargetTrace     = false; // ターゲット内の寄与内訳を毎回printするか
  G4bool   fSaveRndmOnPeak    = false; // ピーク窓に入ったら rndm を保存するか
  G4double fPeakEmin_MeV      = 6.0;   // MeV
  G4double fPeakEmax_MeV      = 10.0;   // MeV
  G4String fRndmDir           = ".";   // 保存先
  G4String fRndmPrefix        = "rndm_evt"; // ファイル名prefix
};

#endif
