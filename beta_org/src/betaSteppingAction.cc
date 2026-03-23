#include "betaSteppingAction.hh"
#include "betaEventAction.hh"

#include "G4Step.hh"
#include "G4VPhysicalVolume.hh"
#include "G4TouchableHandle.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

#include "G4VProcess.hh" // ★これを追加

// #include "betaLineage.hh"

void betaSteppingAction::UserSteppingAction(const G4Step *step)
{
  if (!step)
    return;

  const auto *pre = step->GetPreStepPoint();
  const auto *post = step->GetPostStepPoint();
  const auto *trk = step->GetTrack();
  if (!pre || !post || !trk)
    return;

  const auto *prePV = pre->GetPhysicalVolume();
  const auto *postPV = post->GetPhysicalVolume();
  if (!prePV || !postPV)
    return;

  if (postPV->GetName() == "cellPhysical" && prePV->GetName() != "cellPhysical")
  {

    const auto postTH = post->GetTouchableHandle();
    const int copyNo = postTH ? postTH->GetCopyNumber() : -1;

    const int tid = trk->GetTrackID();
    const int pdg = trk->GetDefinition()->GetPDGEncoding();
    const auto *cp = trk->GetCreatorProcess();
    const auto creator = cp ? cp->GetProcessName() : "Primary";

    // if (trk->GetKineticEnergy() / MeV > 1.0)
    // {
    //   G4cout << "[CAL-ENTER] copyNo=" << copyNo << "\n"
    //         //  << " tid=" << tid
    //          << " pid=" << trk->GetParentID() << "\n"
    //          << " pdg=" << pdg << "\n"
    //          << " Ek=" << trk->GetKineticEnergy() / MeV << " MeV" << "\n"
    //          << " creator=" << creator << "\n"
    //          << " from=" << prePV->GetName() << "\n"
    //          << " to=" << postPV->GetName() << "\n"
    //          << " pos=(" << post->GetPosition().x() / mm << ","
    //          << post->GetPosition().y() / mm << ","
    //          << post->GetPosition().z() / mm << ")mm"
    //          << "\n";
    // }

    // betaLineage::Instance()->DumpLineage(tid, 8);
  }

  const auto edep = step->GetTotalEnergyDeposit();
  if (edep <= 0.)
    return;

  const auto &name = prePV->GetName();
  if (name == "cellPhysical")
  {
    fEvent->AddEdepCell(edep);
  }
  else if (name == "Target")
  {
    fEvent->AddEdepTarget(edep);
  }
}