// calRunAction.cc
#include "calRunAction.hh"
#include "calAnalysis.hh"
#include "calConstant.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"

#include "TROOT.h"
#include "TH1.h"
#include "TTree.h"
#include "TFile.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

calRunAction::calRunAction()
 : G4UserRunAction()
{ 
  // initialize root
  gROOT-> Reset();
  
  // create root file
  hfile = new TFile("test.root","RECREATE", "Geant4 ROOT analysis");

  // define histo
  Cal_dE = new TH1D("h_dE", "Total enegy deposit", 200, 0., 200.);

  // define tree
  tree = new TTree("tree","tree");

  G4int dummy;
  tree->Branch("eveno"		,&dummy,      "eveno/I");

  G4double dummy_beam;
  //tree->Branch("beamx"	,&dummy_beam, "beamx/D");
  //tree->Branch("beamy"	,&dummy_beam, "beamy/D");
  //tree->Branch("beamz"	,&dummy_beam, "beamz/D");
  //tree->Branch("beamMom"	,dummy_beam, "beamMom/D");
  tree->Branch("beamE"		,&dummy_beam, "beamE/D");
  tree->Branch("decayFlag"	,&dummy, "decayFlag/I");


  G4double dummy_data;
  tree->Branch("dEtot_Cal"	,&dummy_data, "dEtot_Cal/D");
  tree->Branch("dEtot_Scinti"	,&dummy_data, "dEtot_Scinti/D");
  tree->Branch("dEtot_Target"	,&dummy_data, "dEtot_Target/D");
  tree->Branch("multiHit"	,&dummy,      "multiHit/I");

  G4int cell_data_int[nb_cryst];
  G4double cell_data[nb_cryst];
  G4String tmp;
  std::ostringstream oss;

  oss << "dE[" << nb_cryst << "]/D";
  tmp=oss.str();
  tree->Branch("dE",cell_data,tmp);
  
  oss.str("");
  oss << "dE_Scinti[" << nb_cryst << "]/D";
  tmp=oss.str();
  tree->Branch("dE_Scinti",cell_data,tmp);

  //oss.str("");
  //oss << "FirstHitPosX[" << nb_cryst << "]/D";
  //tmp=oss.str();
  //tree->Branch("FirstHitPosX",cell_data,tmp);
  //
  //oss.str("");
  //oss << "FirstHitPosY[" << nb_cryst << "]/D";
  //tmp=oss.str();
  //tree->Branch("FirstHitPosY",cell_data,tmp);
  //
  //oss.str("");
  //oss << "FirstHitPosZ[" << nb_cryst << "]/D";
  //tmp=oss.str();
  //tree->Branch("FirstHitPosZ",cell_data,tmp);

  oss.str("");
  oss << "FirstHitTime[" << nb_cryst << "]/D";
  tmp=oss.str();
  tree->Branch("FirstHitTime",cell_data,tmp);

  //oss.str("");
  //oss << "cellNoFlag[" << nb_cryst << "]/I";
  //tmp=oss.str();
  //tree->Branch("cellNoFlag",cell_data_int,tmp);


//  // set printing event number per each event
//  G4RunManager::GetRunManager()->SetPrintProgress(1);     
//
//  // Create analysis manager
//  // The choice of analysis technology is done via selectin of a namespace
//  // in calAnalysis.hh
//  auto analysisManager = G4AnalysisManager::Instance();
//  G4cout << "Using " << analysisManager->GetType() << G4endl;
//
//  // Create directories 
//  //analysisManager->SetHistoDirectoryName("histograms");
//  //analysisManager->SetNtupleDirectoryName("ntuple");
//  analysisManager->SetVerboseLevel(1);
//  analysisManager->SetNtupleMerging(true);
//    // Note: merging ntuples is available only with Root output
//
//  // Book histograms, ntuple
//  //
//  
//  // Creating histograms
//  analysisManager->CreateH1("Eabs","Edep in absorber", 100, 0., 800*MeV);
//  analysisManager->CreateH1("Egap","Edep in gap", 100, 0., 100*MeV);
//  analysisManager->CreateH1("Labs","trackL in absorber", 100, 0., 1*m);
//  analysisManager->CreateH1("Lgap","trackL in gap", 100, 0., 50*cm);
//
//  // Creating ntuple
//  //
//  analysisManager->CreateNtuple("tree", "Edep and TrackL");
//  
//  //analysisManager->CreateNtupleDColumn("Eabs");
//  //analysisManager->CreateNtupleDColumn("Egap");
//
//  std::stringstream ss;
//  for(G4int i= 0; i<nb_cryst; i++){
//    ss.str("");
//    ss << "Eabs" << i;
//    analysisManager->CreateNtupleDColumn(ss.str().c_str());
//  }
//  analysisManager->CreateNtupleDColumn("Labs");
//  analysisManager->CreateNtupleDColumn("Lgap");
//  analysisManager->FinishNtuple();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

calRunAction::~calRunAction()
{
  hfile->Write(); 
  hfile->Close(); 
  delete hfile;
  //delete G4AnalysisManager::Instance();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void calRunAction::BeginOfRunAction(const G4Run* /*run*/)
{ 
  //inform the runManager to save random number seed
  //G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  
  // Get analysis manager
  //auto analysisManager = G4AnalysisManager::Instance();

  // Open an output file
  //G4String fileName = "test";
  //analysisManager->OpenFile(fileName);

  tree-> Reset();

  RunSS_Time = time(NULL);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void calRunAction::EndOfRunAction(const G4Run* /*run*/)
{
  // print histogram statistics
  //
//  auto analysisManager = G4AnalysisManager::Instance();
//  if ( analysisManager->GetH1(1) ) {
//    G4cout << G4endl << " ----> print histograms statistic ";
//    if(isMaster) {
//      G4cout << "for the entire run " << G4endl << G4endl; 
//    }
//    else {
//      G4cout << "for the local thread " << G4endl << G4endl; 
//    }
//    
//    G4cout << " EAbs : mean = " 
//       << G4BestUnit(analysisManager->GetH1(0)->mean(), "Energy") 
//       << " rms = " 
//       << G4BestUnit(analysisManager->GetH1(0)->rms(),  "Energy") << G4endl;
//    
//    G4cout << " EGap : mean = " 
//       << G4BestUnit(analysisManager->GetH1(1)->mean(), "Energy") 
//       << " rms = " 
//       << G4BestUnit(analysisManager->GetH1(1)->rms(),  "Energy") << G4endl;
//    
//    G4cout << " LAbs : mean = " 
//      << G4BestUnit(analysisManager->GetH1(2)->mean(), "Length") 
//      << " rms = " 
//      << G4BestUnit(analysisManager->GetH1(2)->rms(),  "Length") << G4endl;
//
//    G4cout << " LGap : mean = " 
//      << G4BestUnit(analysisManager->GetH1(3)->mean(), "Length") 
//      << " rms = " 
//      << G4BestUnit(analysisManager->GetH1(3)->rms(),  "Length") << G4endl;
//  }

  // save histograms & ntuple
  //
//  analysisManager->Write();
//  analysisManager->CloseFile();

  G4cout << G4endl <<" ** End of Run **" << G4endl;
  RunE_Time = time(NULL) - RunSS_Time;
  G4cout << "Time of This Run = " << RunE_Time << " sec" << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
