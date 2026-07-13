#include "betaRunAction.hh"
#include "BetaConfig.hh"

#include "G4AnalysisManager.hh"
#include "G4Run.hh"

// --- ABLA coupling 用 ---
#include "G4HadronicInteraction.hh"
#include "G4HadronicInteractionRegistry.hh"
#include "G4INCLXXInterface.hh"
#include "G4INCLXXInterfaceStore.hh"
#include "G4AblaInterface.hh"

#include "Constant.hh"
#include "G4SystemOfUnits.hh"
#include "G4Threading.hh"

betaRunAction::betaRunAction(
    std::vector<double> *dEvec, std::vector<int> *pidVec,
    std::vector<double> *thDE, std::vector<double> *thTime,
    std::vector<double> *thPath, std::vector<double> *tlcDE,
    std::vector<double> *tlcCherenkovTime, std::vector<double> *tlcPath,
    std::vector<double> *tlcCherenkovPath,
    std::vector<double> *tlcCherenkovExpectedPhotons)
    : G4UserRunAction(),
      fFileName(BetaConfig::Instance().Output()),
      fDEPtr(dEvec),
      fPIDPtr(pidVec),
      fDummyDE(BetaConfig::Instance().NCells(), 0.0),
      fDummyPID(BetaConfig::Instance().NCells(), 0),
      fTHDEPtr(thDE),
      fTHTimePtr(thTime),
      fTHPathPtr(thPath),
      fTLCDEPtr(tlcDE),
      fTLCCherenkovTimePtr(tlcCherenkovTime),
      fTLCPathPtr(tlcPath),
      fTLCCherenkovPathPtr(tlcCherenkovPath),
      fTLCCherenkovExpectedPhotonsPtr(tlcCherenkovExpectedPhotons),
      fDummyTH(nSegTH, 0.0),
      fDummyTLC(nSegTLC, 0.0)
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
  if (!fTHDEPtr)
    fTHDEPtr = &fDummyTH;
  if (!fTHTimePtr)
    fTHTimePtr = &fDummyTH;
  if (!fTHPathPtr)
    fTHPathPtr = &fDummyTH;
  if (!fTLCDEPtr)
    fTLCDEPtr = &fDummyTLC;
  if (!fTLCCherenkovTimePtr)
    fTLCCherenkovTimePtr = &fDummyTLC;
  if (!fTLCPathPtr)
    fTLCPathPtr = &fDummyTLC;
  if (!fTLCCherenkovPathPtr)
    fTLCCherenkovPathPtr = &fDummyTLC;
  if (!fTLCCherenkovExpectedPhotonsPtr)
    fTLCCherenkovExpectedPhotonsPtr = &fDummyTLC;

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

  // ---- ntuple 3: calarr (geometry-sized) ----
  man->CreateNtuple("calarr", "dE vector + pid vector per event");
  man->CreateNtupleIColumn("eventID");
  man->CreateNtupleDColumn("dE_MeV", *fDEPtr); // vector<double>
  man->CreateNtupleIColumn("pid", *fPIDPtr);   // vector<int>
  man->FinishNtuple();

  // ---- ntuple 4: run metadata ----
  man->CreateNtuple("runmeta", "simulation configuration");
  man->CreateNtupleIColumn("nLayer");
  man->CreateNtupleIColumn("nSector");
  man->CreateNtupleIColumn("segmentationMode");
  man->CreateNtupleIColumn("physicsFlag");
  man->CreateNtupleIColumn("writeCalHit");
  man->CreateNtupleDColumn("thetaMin_deg");
  man->CreateNtupleDColumn("thetaMax_deg");
  man->CreateNtupleDColumn("rMin_cm");
  man->CreateNtupleDColumn("thickness_cm");
  man->CreateNtupleDColumn("neutronScale");
  man->CreateNtupleDColumn("inelasticBias");
  man->CreateNtupleDColumn("pionInelasticXSScale");
  man->CreateNtupleDColumn("seed");
  man->CreateNtupleSColumn("segmentation");
  man->CreateNtupleSColumn("primary");
  man->CreateNtupleSColumn("output");
  man->CreateNtupleIColumn("nSegTH");
  man->CreateNtupleIColumn("nSegTLC");
  man->CreateNtupleDColumn("tlcRefractiveIndex");
  man->CreateNtupleDColumn("tlcLambdaMin_nm");
  man->CreateNtupleDColumn("tlcLambdaMax_nm");
  man->CreateNtupleSColumn("tlcResponseModel");
  man->CreateNtupleDColumn("tlcCollectionEfficiencyApplied");
  man->FinishNtuple();

  // ---- ntuple 5: TH plastic hodoscope event vectors ----
  man->CreateNtuple("th", "TH plastic hodoscope per event");
  man->CreateNtupleIColumn("eventID");
  man->CreateNtupleDColumn("dE_MeV", *fTHDEPtr);
  man->CreateNtupleDColumn("time_ns", *fTHTimePtr);
  man->CreateNtupleDColumn("chargedPath_mm", *fTHPathPtr);
  man->FinishNtuple();

  // ---- ntuple 6: TLC Lucite Cherenkov event vectors ----
  man->CreateNtuple("tlc", "TLC Lucite response per event");
  man->CreateNtupleIColumn("eventID");
  man->CreateNtupleDColumn("dE_truth_MeV", *fTLCDEPtr);
  man->CreateNtupleDColumn("cherenkovTime_ns", *fTLCCherenkovTimePtr);
  man->CreateNtupleDColumn("chargedPath_truth_mm", *fTLCPathPtr);
  man->CreateNtupleDColumn("cherenkovPath_mm", *fTLCCherenkovPathPtr);
  man->CreateNtupleDColumn("cherenkovExpectedPhotons",
                          *fTLCCherenkovExpectedPhotonsPtr);
  man->FinishNtuple();
}

betaRunAction::~betaRunAction()
{
  // G4TaskRunManager clears the thread-local analysis managers at shutdown.
}

void betaRunAction::BeginOfRunAction(const G4Run *)
{
  auto *man = G4AnalysisManager::Instance();
  man->OpenFile(fFileName);

  G4bool fillMetadata = true;
#ifdef G4MULTITHREADED
  fillMetadata = !IsMaster() && G4Threading::G4GetThreadId() == 0;
#endif
  if (fillMetadata)
  {
    const auto &config = BetaConfig::Instance();
    man->FillNtupleIColumn(4, 0, config.NLayer());
    man->FillNtupleIColumn(4, 1, config.NSector());
    man->FillNtupleIColumn(4, 2, config.SegmentationMode());
    man->FillNtupleIColumn(4, 3, PhysicsFlag);
    man->FillNtupleIColumn(4, 4, config.WriteCalHit() ? 1 : 0);
    man->FillNtupleDColumn(4, 5, thetaMin);
    man->FillNtupleDColumn(4, 6, thetaMax);
    man->FillNtupleDColumn(4, 7, Rmin / cm);
    man->FillNtupleDColumn(4, 8, absoThickness / cm);
    man->FillNtupleDColumn(4, 9, NeutronScale);
    man->FillNtupleDColumn(4, 10, InelasticBias);
    man->FillNtupleDColumn(4, 11, PionInelasticXSScale);
    man->FillNtupleDColumn(4, 12, static_cast<double>(config.Seed()));
    man->FillNtupleSColumn(4, 13, config.Segmentation());
    man->FillNtupleSColumn(4, 14, config.Primary());
    man->FillNtupleSColumn(4, 15, config.Output());
    man->FillNtupleIColumn(4, 16, nSegTH);
    man->FillNtupleIColumn(4, 17, nSegTLC);
    man->FillNtupleDColumn(4, 18, TLCRefractiveIndex);
    man->FillNtupleDColumn(4, 19, TLCLambdaMin / nm);
    man->FillNtupleDColumn(4, 20, TLCLambdaMax / nm);
    man->FillNtupleSColumn(4, 21, "analytic_frank_tamm_no_transport");
    man->FillNtupleDColumn(4, 22, 0.0);
    man->AddNtupleRow(4);
  }

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