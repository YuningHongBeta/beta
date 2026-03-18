#include "betaRunAction.hh"

#include "G4AnalysisManager.hh"
#include "G4Run.hh"

// --- ABLA coupling 用 ---
#include "G4HadronicInteraction.hh"
#include "G4HadronicInteractionRegistry.hh"
#include "G4INCLXXInterface.hh"
#include "G4INCLXXInterfaceStore.hh"
#include "G4AblaInterface.hh"

#include "Constant.hh"

betaRunAction::betaRunAction(std::vector<double> *dEvec,
                             std::vector<int> *pidVec)
    : G4UserRunAction(),
      fDEPtr(dEvec),
      fPIDPtr(pidVec),
      fDummyDE(225, 0.0),
      fDummyPID(225, 0)
{
  auto *man = G4AnalysisManager::Instance();
  man->SetVerboseLevel(1);

  // ROOT 出力（Geant4 が ROOT 対応ビルドされている前提）
  man->SetDefaultFileType("root");
  man->SetNtupleMerging(true);

  // nullptr ならダミーに差し替え
  if (!fDEPtr)
    fDEPtr = &fDummyDE;
  if (!fPIDPtr)
    fPIDPtr = &fDummyPID;

  // ---- ntuple 0: evt summary ----
  man->CreateNtuple("evt", "event summary");
  man->CreateNtupleIColumn("eventID");        // col0
  man->CreateNtupleDColumn("EdepCell_MeV");   // col1
  man->CreateNtupleDColumn("EdepTarget_MeV"); // col2
  man->FinishNtuple();

  // ---- ntuple 1: calhit (per segment hit list) ----
  man->CreateNtuple("calhit", "calorimeter hits (per segment per event)");
  man->CreateNtupleIColumn("eventID");
  man->CreateNtupleIColumn("copyNo");
  man->CreateNtupleDColumn("t_ns");
  man->CreateNtupleDColumn("dE_MeV");
  man->CreateNtupleIColumn("pdg");
  man->CreateNtupleDColumn("px_MeV");
  man->CreateNtupleDColumn("py_MeV");
  man->CreateNtupleDColumn("pz_MeV");
  man->CreateNtupleSColumn("creator");
  man->CreateNtupleIColumn("originType"); // 0:unknown, 1:evap, 2:fission
  man->FinishNtuple();

  // ---- ntuple 2: target ----
  man->CreateNtuple("target", "target (per event)");
  man->CreateNtupleIColumn("eventID");
  man->CreateNtupleDColumn("t_ns");
  man->CreateNtupleDColumn("dE_MeV");
  man->CreateNtupleIColumn("pdg");
  man->CreateNtupleDColumn("px_MeV");
  man->CreateNtupleDColumn("py_MeV");
  man->CreateNtupleDColumn("pz_MeV");
  man->FinishNtuple();

  // ---- ntuple 3: calarr (225 fixed) ----
  man->CreateNtuple("calarr", "dE vector(225) + pid vector(225) per event");
  man->CreateNtupleIColumn("eventID");
  man->CreateNtupleDColumn("dE_MeV", *fDEPtr); // vector<double>
  man->CreateNtupleIColumn("pid", *fPIDPtr);   // vector<int>
  man->FinishNtuple();
}

betaRunAction::~betaRunAction()
{
  // G4cout << "DEBUG: ~betaRunAction start" << G4endl;
  // delete G4AnalysisManager::Instance(); 
  // G4cout << "DEBUG: ~betaRunAction end" << G4endl;
}

void betaRunAction::BeginOfRunAction(const G4Run *)
{
  G4AnalysisManager::Instance()->OpenFile(fFileName);

  if (PhysicsFlag == 2)
  {
    // Get hold of pointers to the INCL++ model interface 
    std::vector<G4HadronicInteraction *> interactions = G4HadronicInteractionRegistry::Instance()->FindAllModels(G4INCLXXInterfaceStore::GetInstance()->getINCLXXVersionName());
    for (std::vector<G4HadronicInteraction *>::const_iterator iInter = interactions.begin(), e = interactions.end();
         iInter != e; ++iInter)
    {
      G4INCLXXInterface *theINCLInterface = static_cast<G4INCLXXInterface *>(*iInter);
      if (theINCLInterface)
      {
        // Instantiate the ABLA model
        G4HadronicInteraction *interaction = G4HadronicInteractionRegistry::Instance()->FindModel("ABLA");
        G4AblaInterface *theAblaInterface = static_cast<G4AblaInterface *>(interaction);
        if (!theAblaInterface)
          theAblaInterface = new G4AblaInterface;
        // Couple INCL++ to ABLA
        // G4cout << "Coupling INCLXX to ABLA" << G4endl;
        theINCLInterface->SetDeExcitation(theAblaInterface);
      }
    }
  }
}

void betaRunAction::EndOfRunAction(const G4Run *)
{
  auto *man = G4AnalysisManager::Instance();
  man->Write();
  man->CloseFile();
}