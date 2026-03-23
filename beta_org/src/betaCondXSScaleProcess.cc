#include "betaCondXSScaleProcess.hh"

#include "G4Track.hh"
#include "G4Step.hh"
#include "G4VParticleChange.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"

betaCondXSScaleProcess::betaCondXSScaleProcess(const G4String& name,
                                               G4VProcess* wrapped,
                                               double xsFactor,
                                               std::vector<G4String> targetLVNames,
                                               std::vector<G4String> targetMatNames)
: G4VProcess(name),
  fWrapped(wrapped),
  fXSFactor(xsFactor),
  fTargetLVs(std::move(targetLVNames)),
  fTargetMats(std::move(targetMatNames))
{
}

betaCondXSScaleProcess::~betaCondXSScaleProcess()
{
  delete fWrapped;
  fWrapped = nullptr;
}

G4bool betaCondXSScaleProcess::IsApplicable(const G4ParticleDefinition& p)
{
  return fWrapped ? fWrapped->IsApplicable(p) : false;
}

void betaCondXSScaleProcess::PreparePhysicsTable(const G4ParticleDefinition& p)
{
  if (fWrapped) fWrapped->PreparePhysicsTable(p);
}

void betaCondXSScaleProcess::BuildPhysicsTable(const G4ParticleDefinition& p)
{
  if (fWrapped) fWrapped->BuildPhysicsTable(p);
}

void betaCondXSScaleProcess::StartTracking(G4Track* trk)
{
  if (fWrapped) fWrapped->StartTracking(trk);
}

void betaCondXSScaleProcess::EndTracking()
{
  if (fWrapped) fWrapped->EndTracking();
}

bool betaCondXSScaleProcess::MatchName(const G4String& s, const std::vector<G4String>& list) const
{
  for (const auto& x : list) {
    if (s == x) return true;
  }
  return false;
}

bool betaCondXSScaleProcess::InTarget(const G4Track& track) const
{
  // LV名で限定（logical volume名）
  if (!fTargetLVs.empty()) {
    const auto* pv = track.GetVolume();
    const auto* lv = pv ? pv->GetLogicalVolume() : nullptr;
    if (!lv) return false;
    if (!MatchName(lv->GetName(), fTargetLVs)) return false;
  }

  // 材質名で限定
  if (!fTargetMats.empty()) {
    const auto* mat = track.GetMaterial();
    if (!mat) return false;
    if (!MatchName(mat->GetName(), fTargetMats)) return false;
  }

  // どっちも指定なしは「全空間」になるので注意
  return true;
}

G4double betaCondXSScaleProcess::PostStepGetPhysicalInteractionLength(
    const G4Track& track, G4double previousStepSize, G4ForceCondition* condition)
{
  if (!fWrapped) return DBL_MAX;

  G4double L = fWrapped->PostStepGetPhysicalInteractionLength(track, previousStepSize, condition);

  // XS×f  <=>  mean free path / f
  if (fXSFactor > 0.0 && fXSFactor != 1.0) {
    if (InTarget(track)) {
      L /= fXSFactor;
    }
  }
  return L;
}

G4VParticleChange* betaCondXSScaleProcess::PostStepDoIt(const G4Track& track, const G4Step& step)
{
  return fWrapped ? fWrapped->PostStepDoIt(track, step) : nullptr;
}

G4double betaCondXSScaleProcess::AlongStepGetPhysicalInteractionLength(
    const G4Track& track, G4double previousStepSize, G4double currentMinimumStep,
    G4double& proposedSafety, G4GPILSelection* selection)
{
  return fWrapped
    ? fWrapped->AlongStepGetPhysicalInteractionLength(track, previousStepSize, currentMinimumStep,
                                                     proposedSafety, selection)
    : DBL_MAX;
}

G4VParticleChange* betaCondXSScaleProcess::AlongStepDoIt(const G4Track& track, const G4Step& step)
{
  return fWrapped ? fWrapped->AlongStepDoIt(track, step) : nullptr;
}

G4double betaCondXSScaleProcess::AtRestGetPhysicalInteractionLength(
    const G4Track& track, G4ForceCondition* condition)
{
  return fWrapped ? fWrapped->AtRestGetPhysicalInteractionLength(track, condition) : DBL_MAX;
}

G4VParticleChange* betaCondXSScaleProcess::AtRestDoIt(const G4Track& track, const G4Step& step)
{
  return fWrapped ? fWrapped->AtRestDoIt(track, step) : nullptr;
}