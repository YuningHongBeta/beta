#include "G4RunManagerFactory.hh"
#include "G4MTRunManager.hh"
#include "G4UImanager.hh"
#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"
#include "Randomize.hh"

#include "FTFP_BERT.hh"

#include "betaPhysicsList.hh"
#include "betaDetectorConstruction.hh"
#include "betaActionInitialization.hh"
#include "BetaConfig.hh"

int main(int argc, char **argv)
{
  const auto &config = BetaConfig::Instance();
  G4Random::setTheSeed(config.Seed());

  auto *runManager =
      G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);
#ifdef G4MULTITHREADED
  if (auto *mtRunManager = dynamic_cast<G4MTRunManager *>(runManager))
    mtRunManager->SetNumberOfThreads(config.Threads());
#endif

  G4cout << "[BETA-CONFIG] primary=" << config.Primary()
         << " geometry=" << config.Geometry()
         << " cells=" << config.NLayer() << "x" << config.NSector()
         << " segmentation=" << config.Segmentation()
         << " theta_deg=" << config.ThetaMinDeg() << ".." << config.ThetaMaxDeg()
         << " r_cm=" << config.RMinCm() << ".."
         << config.RMinCm() + config.ThicknessCm()
         << " bgo_z_offset_cm=" << config.BgoZOffsetCm()
         << " photon_counter=" << config.PhotonCounter()
         << " output=" << config.Output()
         << " seed=" << config.Seed()
         << " threads=" << config.Threads() << G4endl;

  runManager->SetUserInitialization(new betaDetectorConstruction());
  if (PhysicsFlag == 1)
  {
    runManager->SetUserInitialization(new FTFP_BERT());
  }
  else
  {
    runManager->SetUserInitialization(new betaPhysicsList());
  }
  runManager->SetUserInitialization(new betaActionInitialization());

  // debug
  // runManager->SetNumberOfThreads(1);
  runManager->Initialize();

  G4VisExecutive *visManager = nullptr;
  if (argc == 1)
  {
    visManager = new G4VisExecutive;
    visManager->Initialize();
  }

  auto *UImanager = G4UImanager::GetUIpointer();

  if (argc == 1)
  {
    auto *ui = new G4UIExecutive(argc, argv);
    auto status = UImanager->ApplyCommand("/control/execute macros/init_vis.mac");
    if (status != 0)
    {
      G4cerr << "Failed to execute init_vis.mac\n";
    }
    ui->SessionStart();
    delete ui;
  }
  else
  {
    G4String cmd = "/control/execute ";
    UImanager->ApplyCommand(cmd + argv[1]);
  }

  delete visManager;
  delete runManager;
  return 0;
}
