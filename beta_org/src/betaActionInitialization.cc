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
    const_cast<std::vector<double>*>(&eventAction->GetDE_MeV()),
    const_cast<std::vector<int>*>(&eventAction->GetPID())
  ));

  SetUserAction(new betaSteppingAction(eventAction));
  SetUserAction(new betaTrackingAction());
  if (PhysicsFlag >= 4)
  {
    SetUserAction(new betaStackingAction(NeutronScale));
  }
}