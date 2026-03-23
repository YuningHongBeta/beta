#pragma once
#include "G4VUserActionInitialization.hh"

class betaActionInitialization : public G4VUserActionInitialization {
public:
  betaActionInitialization() = default;
  ~betaActionInitialization() override = default;

  void BuildForMaster() const override;
  void Build() const override;
};