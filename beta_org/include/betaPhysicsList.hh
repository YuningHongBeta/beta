#pragma once

#include "G4VModularPhysicsList.hh"

class betaPhysicsList : public G4VModularPhysicsList {
public:
  betaPhysicsList();
  ~betaPhysicsList() override = default;

  void SetCuts() override;
  void ConstructProcess() override;
};