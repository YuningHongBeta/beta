#pragma once

#include "G4VProcess.hh"
#include "G4String.hh"

#include <vector>

class betaCondXSScaleProcess : public G4VProcess
{
public:
  betaCondXSScaleProcess(const G4String& name,
                         G4VProcess* wrapped,
                         double xsFactor,
                         std::vector<G4String> targetLVNames = {},
                         std::vector<G4String> targetMatNames = {});
  ~betaCondXSScaleProcess() override;

  // delegate
  G4bool IsApplicable(const G4ParticleDefinition& p) override;
  void   PreparePhysicsTable(const G4ParticleDefinition& p) override;
  void   BuildPhysicsTable(const G4ParticleDefinition& p) override;
  void   StartTracking(G4Track* trk) override;
  void   EndTracking() override;

  // step methods (discrete)
  G4double PostStepGetPhysicalInteractionLength(
      const G4Track& track,
      G4double previousStepSize,
      G4ForceCondition* condition) override;

  G4VParticleChange* PostStepDoIt(const G4Track& track, const G4Step& step) override;

  // pass-through (not used typically)
  G4double AlongStepGetPhysicalInteractionLength(
      const G4Track& track, G4double previousStepSize, G4double currentMinimumStep,
      G4double& proposedSafety, G4GPILSelection* selection) override;

  G4VParticleChange* AlongStepDoIt(const G4Track& track, const G4Step& step) override;

  G4double AtRestGetPhysicalInteractionLength(
      const G4Track& track, G4ForceCondition* condition) override;

  G4VParticleChange* AtRestDoIt(const G4Track& track, const G4Step& step) override;

private:
  bool InTarget(const G4Track& track) const;
  bool MatchName(const G4String& s, const std::vector<G4String>& list) const;

private:
  G4VProcess* fWrapped = nullptr;  // owned
  double      fXSFactor = 1.0;     // XS倍率 (>1 で増える)
  std::vector<G4String> fTargetLVs;
  std::vector<G4String> fTargetMats;
};