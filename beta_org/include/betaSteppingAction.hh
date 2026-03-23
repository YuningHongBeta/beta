// betaSteppingAction.hh
#pragma once
#include "G4UserSteppingAction.hh"
#include <unordered_set>

class betaEventAction;
class G4Step;

class betaSteppingAction : public G4UserSteppingAction {
public:
  explicit betaSteppingAction(betaEventAction* evt) : fEvent(evt) {}
  ~betaSteppingAction() override = default;

  void UserSteppingAction(const G4Step*) override;

  // eventごとにリセットしたいなら EventAction から呼ぶ
  void PrepareNewEvent() { fEnteredOnce.clear(); }

  // ログON/OFF（必要なら）
  void SetVerbose(bool v) { fVerbose = v; }

private:
  betaEventAction* fEvent = nullptr;

  bool fVerbose = true;

  // tid×copyNo の「初回侵入」だけ出すため
  std::unordered_set<unsigned long long> fEnteredOnce;
};