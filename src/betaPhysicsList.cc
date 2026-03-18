#include "betaPhysicsList.hh"

#include <iomanip>
#include <CLHEP/Units/SystemOfUnits.h>

#include <globals.hh>
#include <G4ios.hh>

#include <G4DecayPhysics.hh>
#include <G4EmStandardPhysics.hh>
#include <G4EmExtraPhysics.hh>
#include <G4IonPhysics.hh>
#include <G4StoppingPhysics.hh>
#include <G4HadronElasticPhysics.hh>
#include <G4NeutronTrackingCut.hh>
#include <G4HadronPhysicsQGSP_BERT.hh>
#include <G4HadronPhysicsFTFP_BERT_HP.hh>
#include <G4HadronPhysicsINCLXX.hh>
#include <G4HadronElasticPhysicsHP.hh>
#include <G4IonINCLXXPhysics.hh>
#include <G4VPhysicsConstructor.hh>
// #include <G4HadronPhysicsShielding.hh>
// #include <G4OpticalPhysics.hh>
#include <G4PionMinus.hh>
#include <G4ProcessManager.hh>
#include <G4HadronicProcess.hh>
#include <G4GenericBiasingPhysics.hh>

// for GEM
#include "G4HadronicInteractionRegistry.hh"
#include "G4INCLXXInterfaceStore.hh"
#include "G4INCLXXInterface.hh"
#include "G4PreCompoundModel.hh"
#include "G4ExcitationHandler.hh"
#include "G4Evaporation.hh"
#include "G4NuclearLevelData.hh"
#include "G4DeexPrecoParameters.hh"

#include <G4HadronPhysicsFTFP_BERT_HP.hh>
#include <G4HadronPhysicsINCLXX.hh>
#include <G4HadronElasticPhysicsHP.hh>
#include <G4IonINCLXXPhysics.hh>
#include <G4VPhysicsConstructor.hh>
#include "G4HadronStoppingProcess.hh"
// #include <G4HadronPhysicsShielding.hh>
// #include <G4OpticalPhysics.hh>
#include <G4PionMinus.hh>
#include <G4ProcessManager.hh>
#include <G4VProcess.hh>
#include <G4HadronicProcess.hh>
#include <G4GenericBiasingPhysics.hh>
#include <G4HadronicInteractionRegistry.hh>
#include <G4INCLXXInterface.hh>
#include <G4ExcitationHandler.hh>
#include <G4Evaporation.hh>
#include <G4PreCompoundModel.hh>
#include <G4INCLXXInterfaceStore.hh>

#include "betaCondXSScaleProcess.hh"   // または .h

#include "Constant.hh"

betaPhysicsList::betaPhysicsList()
{
  // --- ここで de-excitation チャネルを選ぶ（PreInit中に実行される） ---
  auto* deex = G4NuclearLevelData::GetInstance()->GetParameters();
  deex->SetDeexChannelsType(fGEM);      // GEM を使う
  // deex->SetDeexChannelsType(fCombined); // Evaporation+GEM にしたいならこちら

  // verbosity (0: silent, 1: basic, 2+: chatty)
  SetVerboseLevel(0);

  // Default cut value (length). You can tune later.
  defaultCutValue = 0.7 * CLHEP::mm;

  // if (PhysicsFlag == 4 || PhysicsFlag == 5)
  // {
  //   auto *biasingPhysics = new G4GenericBiasingPhysics();
  //   biasingPhysics->Bias("pi-");
  //   RegisterPhysics(biasingPhysics);
  // }

  // --- Electromagnetic ---
  RegisterPhysics(new G4EmStandardPhysics());

  // Extra EM (gamma-nuclear, etc. useful depending on your case)
  RegisterPhysics(new G4EmExtraPhysics());

  // --- Decays ---
  RegisterPhysics(new G4DecayPhysics());

  // --- Hadronics ---
  // RegisterPhysics(new G4HadronElasticPhysics());
  RegisterPhysics(new G4HadronElasticPhysicsHP());
  // RegisterPhysics(new G4HadronPhysicsFTFP_BERT());
  RegisterPhysics(new G4HadronPhysicsINCLXX("hInelasticINCLXX", true, true));

  // --- Stopping (mu-, hadron stopping) ---
  RegisterPhysics(new G4StoppingPhysics());

  // --- Ions ---
  // RegisterPhysics(new G4IonPhysics());
  RegisterPhysics(new G4IonINCLXXPhysics());

  // --- Neutron tracking cut (performance helper) ---
  // RegisterPhysics(new G4NeutronTrackingCut());
}

void betaPhysicsList::SetCuts()
{
  // Use defaultCutValue for gamma, e-, e+ by default
  SetCutsWithDefault();

  // If you want per-particle cuts, you can do e.g.
  SetCutValue(0.1 * CLHEP::mm, "gamma");
  SetCutValue(0.1 * CLHEP::mm, "e-");
  SetCutValue(0.1 * CLHEP::mm, "e+");

  if (verboseLevel > 0)
  {
    DumpCutValuesTable();
  }
}

void betaPhysicsList::ConstructProcess()
{
  // Let the base class construct processes
  G4VModularPhysicsList::ConstructProcess();

  if (PhysicsFlag == 3 || PhysicsFlag == 4 || PhysicsFlag == 5)
  {

    // Appling GEM
    // 2. INCLXXモデルを探す
    auto *reg = G4HadronicInteractionRegistry::Instance();

    // INCL++ が registry に登録している “本当の名前”
    const auto &inclName =
        G4INCLXXInterfaceStore::GetInstance()->getINCLXXVersionName(); // "INCL++ x.y.z" 等  [oai_citation:2‡nuclear.korea.ac.kr](https://nuclear.korea.ac.kr/~ejungwoo/geant4/G4INCLXXInterfaceStore_8cc_source.html)

    // 同名モデルが複数あるのでまとめて取る  [oai_citation:3‡nuclear.korea.ac.kr](https://nuclear.korea.ac.kr/~ejungwoo/geant4/classG4HadronicInteractionRegistry.html)
    auto models = reg->FindAllModels(inclName);

    if (models.empty())
    {
      G4cout << "### Warning: INCLXX model not found in registry. name=\""
             << inclName << "\"\n";
      return;
    }

    auto *handler = new G4ExcitationHandler();
    auto *evaporation = new G4Evaporation();
    evaporation->SetGEMChannel();
    handler->SetEvaporation(evaporation);

    auto *preCompound = new G4PreCompoundModel(handler);

    int nApplied = 0;
    for (auto *m : models)
    {
      if (auto *incl = dynamic_cast<G4INCLXXInterface *>(m))
      {
        incl->SetDeExcitation(preCompound);
        ++nApplied;
      }
    }

    G4cout << "### GEM applied to INCLXX models: " << nApplied << G4endl;

    // // pi- のプロセスマネージャーを取得
    // G4ProcessManager *pManager = G4PionMinus::PionMinus()->GetProcessManager();
    // G4ProcessVector *processList = pManager->GetProcessList();

    // // プロセスリストを走査して非弾性散乱を探す
    // for (size_t i = 0; i < processList->size(); i++)
    // {
    //   G4VProcess *proc = (*processList)[i];

    //   // プロセス名で判定（通常は "pi-Inelastic"）
    //   if (proc->GetProcessName() == "pi-Inelastic")
    //   {

    //     // G4HadronicProcess にキャスト
    //     G4HadronicProcess *hadProc = dynamic_cast<G4HadronicProcess *>(proc);

    //     if (hadProc)
    //     {
    //       // ここで倍率を指定（例：断面積を1.5倍にする）
    //       // 1.8だと多すぎるので1.5に変更
    //       // 1.5だと不足気味なので1.7に変更
    //       // 1.7だと多いっぽいので1.6で試す
    //       G4double scaleFactor = 1.6;
    //       hadProc->MultiplyCrossSectionBy(scaleFactor);

    //       G4cout << "!!! WARNING: pi- Inelastic XS is scaled by "
    //              << scaleFactor << " !!!" << G4endl;
    //     }
    //   }
    // }

    // -----------------------------
    // pi- inelastic XS scaling (conditional)
    // -----------------------------
    auto *pm = G4PionMinus::PionMinus()->GetProcessManager();
    auto *pv = pm->GetProcessList();

    // G4INCLXXInterface * inclxxModel = new G4INCLXXInterface();
    // G4HadronStoppingProcess* piMinusAbsorption = new G4HadronStoppingProcess("PionMinusAbsorption");
    // piMinusAbsorption->RegisterMe(inclxxModel);
    // pm->AddRestProcess(piMinusAbsorption);

    G4VProcess *old = nullptr;

    for (int i = 0; i < (int)pv->size(); ++i)
    {
      auto *p = (*pv)[i];
      if (!p)
        continue;
      if (p->GetProcessName() == "pi-Inelastic")
      {
        old = p;
        break;
      }
    }

    if (!old)
    {
      G4cout << "### WARNING: pi-Inelastic not found; no XS scaling applied.\n";
      return;
    }

    // ProcessManager から外す（old は wrapper が所有して delete する）
    pm->RemoveProcess(old);

    const double scaleFactor = 1.65;

    // ---- 条件をここで指定 ----
    // LV名（logical volume名）で限定：あなたのコードなら TargetLV の名前は "Target"
    std::vector<G4String> targetLVs = {};

    // 材質名で限定：あなたの enriched Li は "Li6_90pct"
    std::vector<G4String> targetMats = {"G4_BGO", "BGO"};

    // LVだけにしたいなら targetMats を {} にする
    // 材質だけにしたいなら targetLVs を {} にする

    auto *scaled =
        new betaCondXSScaleProcess("pi-Inelastic_scaledInBGO",
                                   old,
                                   scaleFactor,
                                   targetLVs,
                                   targetMats);

    pm->AddDiscreteProcess(scaled);

    G4cout << "!!! pi- Inelastic XS scaled by " << scaleFactor
           << " only when (LV in targetLVs) AND (Material in targetMats) !!!\n";
  }
}