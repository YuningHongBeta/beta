#include "HodoscopeSD.hh"

#include "BetaConfig.hh"
#include "Constant.hh"
#include "DetectorResponse.hh"
#include "G4HCofThisEvent.hh"
#include "G4PhysicalConstants.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4TouchableHandle.hh"
#include "G4Track.hh"

#include <algorithm>
#include <cmath>

static_assert(DetectorResponse::THBarZMax - DetectorResponse::THBarZMin ==
                  TLength,
              "TH timing endpoints must match the TH geometry length");

HodoscopeSD::HodoscopeSD(const G4String &name,
                         const G4String &hitsCollectionName,
                         G4bool calculateCherenkov)
    : G4VSensitiveDetector(name),
      fCalculateCherenkov(calculateCherenkov)
{
  collectionName.insert(hitsCollectionName);
}

void HodoscopeSD::Initialize(G4HCofThisEvent *hce)
{
  fHitsCollection =
      new HodoscopeHitsCollection(SensitiveDetectorName, collectionName[0]);
  if (fHCID < 0)
    fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection);
  hce->AddHitsCollection(fHCID, fHitsCollection);
  fHitByCopyNo.clear();
}

G4bool HodoscopeSD::ProcessHits(G4Step *step, G4TouchableHistory *)
{
  if (!step || !step->GetTrack() || !step->GetPreStepPoint())
    return false;

  const auto *track = step->GetTrack();
  const auto *pre = step->GetPreStepPoint();
  const auto *post = step->GetPostStepPoint();
  const G4double stepTime = post
                                ? 0.5 * (pre->GetGlobalTime() +
                                         post->GetGlobalTime())
                                : pre->GetGlobalTime();
  const G4double gateWidth = BetaConfig::Instance().HodoGateWidthNs() * ns;
  if (gateWidth > 0.0 && std::abs(stepTime) > 0.5 * gateWidth)
    return false;
  const auto touchable = pre->GetTouchableHandle();
  if (!touchable)
    return false;
  const G4int copyNo = touchable->GetCopyNumber(0);
  if (copyNo < 0)
    return false;

  const G4double edep = step->GetTotalEnergyDeposit();
  const G4double charge = track->GetDynamicParticle()->GetCharge() / eplus;
  const G4double chargedPath =
      (charge != 0.0) ? step->GetStepLength() : 0.0;

  G4double cherenkovPath = 0.0;
  G4double cherenkovExpectedPhotons = 0.0;
  if (fCalculateCherenkov && chargedPath > 0.0)
  {
    const G4double betaPre = pre->GetBeta();
    const G4double betaPost = post ? post->GetBeta() : betaPre;
    const G4double betaThreshold = 1.0 / TLCRefractiveIndex;
    const auto frankTammFactor = [](G4double beta)
    {
      const G4double betaN = beta * TLCRefractiveIndex;
      return betaN > 1.0 ? 1.0 - 1.0 / (betaN * betaN) : 0.0;
    };

    const G4double factorPre = frankTammFactor(betaPre);
    const G4double factorPost = frankTammFactor(betaPost);
    const G4double factorMid = frankTammFactor(0.5 * (betaPre + betaPost));

    if (factorPre > 0.0 && factorPost > 0.0)
      cherenkovPath = chargedPath;
    else if (factorPre > 0.0 || factorPost > 0.0)
    {
      const G4double betaAbove = factorPre > 0.0 ? betaPre : betaPost;
      const G4double betaBelow = factorPre > 0.0 ? betaPost : betaPre;
      const G4double fraction =
          (betaAbove - betaThreshold) / (betaAbove - betaBelow);
      cherenkovPath = chargedPath * std::clamp(fraction, 0.0, 1.0);
    }

    const G4double charge2 = charge * charge;
    const G4double spectralIntegral =
        1.0 / TLCLambdaMin - 1.0 / TLCLambdaMax;
    cherenkovExpectedPhotons =
        2.0 * pi * fine_structure_const * charge2 *
        spectralIntegral * chargedPath *
        (factorPre + 4.0 * factorMid + factorPost) / 6.0;
  }

  if (edep <= 0.0 && chargedPath <= 0.0 && cherenkovExpectedPhotons <= 0.0)
    return false;

  if (static_cast<size_t>(copyNo) >= fHitByCopyNo.size())
    fHitByCopyNo.resize(copyNo + 1, nullptr);
  auto *&hit = fHitByCopyNo[copyNo];
  if (!hit)
  {
    hit = new HodoscopeHit();
    hit->SetCopyNo(copyNo);
    fHitsCollection->insert(hit);
  }

  if (edep > 0.0)
  {
    const G4double depositTime = post
                                     ? 0.5 * (pre->GetGlobalTime() +
                                              post->GetGlobalTime())
                                     : pre->GetGlobalTime();
    hit->UpdateDepositTime(
        depositTime, track->GetDefinition()->GetPDGEncoding());

    if (!fCalculateCherenkov)
    {
      const G4double depositZ = post
                                    ? 0.5 * (pre->GetPosition().z() +
                                             post->GetPosition().z())
                                    : pre->GetPosition().z();
      const G4double leftDistance = std::clamp(
          depositZ - DetectorResponse::THBarZMin,
          0.0, DetectorResponse::THBarZMax - DetectorResponse::THBarZMin);
      const G4double rightDistance = std::clamp(
          DetectorResponse::THBarZMax - depositZ,
          0.0, DetectorResponse::THBarZMax - DetectorResponse::THBarZMin);
      hit->UpdateTHArrivalTimes(
          depositTime + leftDistance / DetectorResponse::THEffectiveLightSpeed,
          depositTime + rightDistance / DetectorResponse::THEffectiveLightSpeed);
    }
  }

  hit->AddEdep(edep);
  hit->AddChargedPath(chargedPath);
  hit->AddCherenkovPath(cherenkovPath);
  if (cherenkovExpectedPhotons > 0.0)
    hit->UpdateCherenkovTime(pre->GetGlobalTime());
  hit->AddCherenkovExpectedPhotons(cherenkovExpectedPhotons);
  return true;
}
