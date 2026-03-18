#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"

#include "FTFP_BERT.hh"

#include "betaPhysicsList.hh"
#include "betaDetectorConstruction.hh"
#include "betaActionInitialization.hh"

int main(int argc, char **argv)
{
  auto *runManager =
      G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);

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

  auto *visManager = new G4VisExecutive;
  visManager->Initialize();

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