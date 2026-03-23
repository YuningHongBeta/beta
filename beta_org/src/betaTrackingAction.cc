#include "betaTrackingAction.hh"

// #include "betaLineage.hh"
#include "G4Track.hh"

void betaTrackingAction::PreUserTrackingAction(const G4Track* trk)
{
  // betaLineage::Instance()->Record(trk);
}