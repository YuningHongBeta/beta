#ifndef betaTrackingAction_hh
#define betaTrackingAction_hh 1

#include "G4UserTrackingAction.hh"

class G4Track;

class betaTrackingAction : public G4UserTrackingAction {
public:
  betaTrackingAction() = default;
  ~betaTrackingAction() override = default;

  void PreUserTrackingAction(const G4Track* trk) override;
};

#endif