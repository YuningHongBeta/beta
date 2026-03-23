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
                         std::vector<int>*    pidVec = nullptr);

  ~betaRunAction() override;

  void BeginOfRunAction(const G4Run*) override;
  void EndOfRunAction(const G4Run*) override;

private:
  G4String fFileName = "output/beta"; // -> output/beta.root

  std::vector<double>* fDEPtr  = nullptr;
  std::vector<int>*    fPIDPtr = nullptr;

  // nullptr の場合の受け皿
  std::vector<double> fDummyDE;
  std::vector<int>    fDummyPID;
};

#endif