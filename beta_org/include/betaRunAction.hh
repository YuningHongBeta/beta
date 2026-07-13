#ifndef betaRunAction_hh
#define betaRunAction_hh 1

#include "G4UserRunAction.hh"
#include "G4String.hh"
#include <vector>

class G4Run;

class betaRunAction : public G4UserRunAction
{
public:
  // dEvec/pidVec を渡すと、その vector を calarr の branch に直結する
  // (nullptr の場合は内部ダミー vector を使う)
  explicit betaRunAction(std::vector<double>* dEvec = nullptr,
                         std::vector<int>* pidVec = nullptr,
                         std::vector<double>* thDE = nullptr,
                         std::vector<double>* thTime = nullptr,
                         std::vector<double>* thTimeLeft = nullptr,
                         std::vector<double>* thTimeRight = nullptr,
                         std::vector<double>* thTimeLeftMinusRight = nullptr,
                         std::vector<double>* thZReco = nullptr,
                         std::vector<double>* thPathTruth = nullptr,
                         std::vector<double>* tlcDE = nullptr,
                         std::vector<double>* tlcCherenkovTime = nullptr,
                         std::vector<double>* tlcPath = nullptr,
                         std::vector<double>* tlcCherenkovPath = nullptr,
                         std::vector<double>* tlcCherenkovExpectedPhotons = nullptr);

  ~betaRunAction() override;

  void BeginOfRunAction(const G4Run*) override;
  void EndOfRunAction(const G4Run*) override;

private:
  G4String fFileName;

  std::vector<double>* fDEPtr  = nullptr;
  std::vector<int>*    fPIDPtr = nullptr;

  // nullptr の場合の受け皿
  std::vector<double> fDummyDE;
  std::vector<int>    fDummyPID;
  std::vector<double>* fTHDEPtr = nullptr;
  std::vector<double>* fTHTimePtr = nullptr;
  std::vector<double>* fTHTimeLeftPtr = nullptr;
  std::vector<double>* fTHTimeRightPtr = nullptr;
  std::vector<double>* fTHTimeLeftMinusRightPtr = nullptr;
  std::vector<double>* fTHZRecoPtr = nullptr;
  std::vector<double>* fTHPathTruthPtr = nullptr;
  std::vector<double>* fTLCDEPtr = nullptr;
  std::vector<double>* fTLCCherenkovTimePtr = nullptr;
  std::vector<double>* fTLCPathPtr = nullptr;
  std::vector<double>* fTLCCherenkovPathPtr = nullptr;
  std::vector<double>* fTLCCherenkovExpectedPhotonsPtr = nullptr;
  std::vector<double> fDummyTH;
  std::vector<double> fDummyTLC;
};

#endif