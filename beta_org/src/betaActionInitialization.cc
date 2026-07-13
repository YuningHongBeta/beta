#include "betaActionInitialization.hh"

#include "betaPrimaryGeneratorAction.hh"
#include "betaRunAction.hh"
#include "betaEventAction.hh"
#include "betaSteppingAction.hh"
#include "betaTrackingAction.hh"
#include "betaStackingAction.hh"

#include "Constant.hh"

void betaActionInitialization::BuildForMaster() const
{
  // master thread: event は回さないので vector を渡さない
  SetUserAction(new betaRunAction());
}

void betaActionInitialization::Build() const
{
  SetUserAction(new betaPrimaryGeneratorAction());

  auto* eventAction = new betaEventAction();
  SetUserAction(eventAction);

  // eventAction が持つ vector を RunAction の branch に直結
  SetUserAction(new betaRunAction(
    &eventAction->GetDE_MeV(),
    &eventAction->GetPID(),
    &eventAction->GetTHDE_MeV(),
    &eventAction->GetTHTime_ns(),
    &eventAction->GetTHTimeLeft_ns(),
    &eventAction->GetTHTimeRight_ns(),
    &eventAction->GetTHTimeLeftMinusRight_ns(),
    &eventAction->GetTHZReco_mm(),
    &eventAction->GetTHPathTruth_mm(),
    &eventAction->GetTLCDE_MeV(),
    &eventAction->GetTLCCherenkovTime_ns(),
    &eventAction->GetTLCPath_mm(),
    &eventAction->GetTLCCherenkovPath_mm(),
    &eventAction->GetTLCCherenkovExpectedPhotons()
  ));

  SetUserAction(new betaSteppingAction(eventAction));
  SetUserAction(new betaTrackingAction());
  if (PhysicsFlag >= 4)
  {
    SetUserAction(new betaStackingAction(NeutronScale));
  }
}