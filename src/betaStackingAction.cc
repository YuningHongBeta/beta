#include "betaStackingAction.hh"

#include "G4StackManager.hh"
#include "G4Track.hh"
#include "G4Neutron.hh"
#include "G4Gamma.hh"
#include "G4VProcess.hh"
#include "G4ios.hh"

#include "Randomize.hh"
// #include "betaLineage.hh"

betaStackingAction::~betaStackingAction() {
    // G4cout << "DEBUG: ~betaStackingAction" << G4endl;
}

void betaStackingAction::PrepareNewEvent()
{
  fromPiMinus_.clear();
  alreadySplit_.clear();
}

// ★核っぽいγ判定（必要ならログで増減）
static bool IsNuclearGammaLike(const G4Track *t)
{
  const auto *cp = t->GetCreatorProcess();
  if (!cp)
    return false;

  const auto &pn = cp->GetProcessName();

  if (pn == "photonNuclear")
    return true;
  if (pn == "nCapture")
    return true;
  if (pn == "neutronInelastic")
    return true;
  if (pn == "hadInelastic")
    return true;
  if (pn == "RadioactiveDecay")
    return true;

  if (pn.find("Deexc") != std::string::npos)
    return true;
  if (pn.find("deexc") != std::string::npos)
    return true;

  if (pn == "hBertiniCaptureAtRest")
    return true;
  if (pn == "PionMinusAbsorption")
    return true;

  return false;
}

bool betaStackingAction::IsPiMinusSeedProcess(const G4Track *t) const
{
  const auto *cp = t->GetCreatorProcess();
  if (!cp)
    return false;

  const auto &pn = cp->GetProcessName();

  if (pn == "pi-Inelastic")
    return true;
  if (pn == "PionMinusInelastic")
    return true;

  if (pn == "PiMinusAbsorptionAtRest")
    return true;
  if (pn == "PionMinusAbsorptionAtRest")
    return true;

  if (pn == "hBertiniCaptureAtRest")
    return true;
  if (pn == "PionMinusAbsorption")
    return true;
  if (pn == "pi-Inelastic_scaledInBGO")
    return true;

  return false;
}

G4ClassificationOfNewTrack
betaStackingAction::ClassifyNewTrack(const G4Track *t)
{
  // betaLineage::Instance()->Record(t);
  const int tid = t->GetTrackID();
  const int pid = t->GetParentID();

  // クローン判定（無限増殖防止）
  bool isClone = false;
  if (auto *info = dynamic_cast<betaTrackInfo *>(t->GetUserInformation()))
    isClone = info->IsClone();

  // ---- π-起源判定（親フラグ伝搬 + 起点判定）
  bool fromPi = false;

  if (pid > 0)
  {
    auto it = fromPiMinus_.find(pid);
    if (it != fromPiMinus_.end() && it->second)
      fromPi = true;
  }

  if (!fromPi && IsPiMinusSeedProcess(t))
    fromPi = true;

  if (fromPi)
    fromPiMinus_[tid] = true;

  // π-起源以外は何もしない
  if (!fromPi)
    return fUrgent;

  const bool isNeutron = (t->GetDefinition() == G4Neutron::Neutron());
  const bool isGamma = (t->GetDefinition() == G4Gamma::Gamma());

  // -----------------------------
  // ★ gamma：核っぽいγだけ 50% kill
  // -----------------------------
  if (isGamma && IsNuclearGammaLike(t))
  {
    if (G4UniformRand() < 0.5)
      return fKill;
    return fUrgent;
  }

  // -----------------------------
  // ★低E中性子：2回に1回 kill（増やさない）
  // -----------------------------
  if (isNeutron && neutronCloneEkinMin_ > 0.0)
  {
    const double ekin = t->GetKineticEnergy();
    if (ekin < neutronCloneEkinMin_)
    {
      if (G4UniformRand() < 0.5)
        return fKill;
      return fUrgent;
    }
  }

  // ここから下は「スプリット」の話なので clone は除外
  if (isClone)
    return fUrgent;

  // Determine origin type (Evap vs Fission)
  betaTrackInfo::OriginType origin = betaTrackInfo::kUnknown;
  
  // Inherit from parent if available
  if (auto* parentInfo = dynamic_cast<betaTrackInfo*>(t->GetUserInformation())) {
    origin = parentInfo->GetOrigin();
  }

  // If currently unknown, try to determine from creator model
  if (origin == betaTrackInfo::kUnknown) {
    const G4String modelName = t->GetCreatorModelName();
    if (modelName.find("Evaporation") != std::string::npos) {
      origin = betaTrackInfo::kEvaporation;
    } else if (modelName.find("Fission") != std::string::npos || modelName.find("fission") != std::string::npos) {
      origin = betaTrackInfo::kFission;
    }
  }

  // --- スプリット対象は中性子のみ（今の方針を維持）
  if (!isNeutron) {
      // Even if not splitting, we want to attach info if it's missing or needs update?
      // For now, only relevant if we are tracking this info for analysis.
      // But standard tracks might not have betaTrackInfo unless we attach it.
      // If t doesn't have UserInfo, we can attach one. 
      if (!t->GetUserInformation()) {
          // const_cast is nasty but necessary if we want to attach info in StackingAction for the primary track object
          // actually ClassifyNewTrack gets a const G4Track*, so we cannot modify it easily without const_cast.
          // However, standard G4 practice in StackingAction allows setting user info on the track if it's in the stack.
          // Wait, we can't modify 't' here safely if it's const. 
          // But we are returning a Classification.
          
          // Actually, for the *clones*, we create new tracks, so we can set info there.
          // For the original track 't', if we want to tag same-track, we might need to modify it.
          // Let's assume we care mostly about the clones or we rely on the fact that we might not be able to tag the *original* successfully here without const_cast.
          // Usage of const_cast is typical in StackingAction for this purpose.
          G4Track* nonConstT = const_cast<G4Track*>(t);
          nonConstT->SetUserInformation(new betaTrackInfo(fromPi, isClone, origin));
      }
      return fUrgent;
  }

  // 同一trackを二重にスプリットしない
  if (!alreadySplit_.insert(tid).second)
    return fUrgent;
  // if (G4UniformRand()<1)
    // return fKill; // temp

  // If original track doesn't have info, attach it now
  if (!t->GetUserInformation()) {
      G4Track* nonConstT = const_cast<G4Track*>(t);
      nonConstT->SetUserInformation(new betaTrackInfo(fromPi, isClone, origin));
  }

  const G4double w0 = t->GetWeight();

  for (int i = 0; i < nClones_; ++i)
  {
    auto *c = new G4Track(*t);
    c->SetParentID(t->GetParentID());
    c->SetWeight(w0); // ★安全のため明示
    c->SetUserInformation(new betaTrackInfo(/*fromPiMinus=*/true, /*isClone=*/true, origin));
    stackManager->PushOneTrack(c);
  }

  return fUrgent;
}