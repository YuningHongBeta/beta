// calEventAction.cc
#include "calEventAction.hh"
//#include "ParamMan.hh"
#include "calRunAction.hh"
#include "calAnalysis.hh"
#include "calConstant.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>

#include "TROOT.h"
#include "TTree.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

calEventAction::calEventAction()
 : G4UserEventAction(),
   evID(0),
   beamx(0.),beamy(0.),beamz(0.),
   beamMom(0.),beamE(0.),
   //paramMan(0),
   fEnergyAbs(0.),
   fTrackLAbs(0.)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

calEventAction::~calEventAction()
{}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void calEventAction::BeginOfEventAction(const G4Event* /*event*/)
{  
  // initialisation per event
  evID = 0;
  beamx = 0.;
  beamy = 0.;
  beamz = 0.;
  beamMom = 0.;
  beamE = 0.;
  fEnergyAbs = 0.;
  fTrackLAbs = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void calEventAction::EndOfEventAction(const G4Event* event)
{
  TTree* tree = (TTree*)gROOT-> FindObject("tree");
  tree->SetBranchAddress("eveno"  ,&evID);

  evID = event->GetEventID();
  //analysisManager->FillNtupleDColumn(0, fEnergyAbs);
  tree->Fill();
  
  
  if( evID%10000 == 0 ){
    G4cout << "End of Event# : " << evID << G4endl;
    //getchar();
    }
  
}
//  // Accumulate statistics
//  //
//
//  // get analysis manager
//  auto analysisManager = G4AnalysisManager::Instance();
//
//  // fill histograms
//  analysisManager->FillH1(0, fEnergyAbs);
//  analysisManager->FillH1(1, fEnergyGap);
//  analysisManager->FillH1(2, fTrackLAbs);
//  analysisManager->FillH1(3, fTrackLGap);
//
//  
//  // fill ntuple
//  analysisManager->FillNtupleDColumn(0, fEnergyAbs);
//  analysisManager->FillNtupleDColumn(1, fEnergyGap);
//  analysisManager->FillNtupleDColumn(2, fTrackLAbs);
//  analysisManager->FillNtupleDColumn(3, fTrackLGap);
//  for(G4int i=0; i<nb_cryst; i++){
//    analysisManager->FillNtupleDColumn(i, fEnergyAbsbyCell[i]);
//    //G4cout << fEnergyAbsbyCell[i] << G4endl;
//    //getchar();
//  }
//
//  analysisManager->AddNtupleRow();  
  
  // Print per event (modulo n)
  //
//  auto eventID = event->GetEventID();
//  auto printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
//  if ( ( printModulo > 0 ) && ( eventID % printModulo == 0 ) ) {
//    G4cout << "---> End of event: " << eventID << G4endl;     
//
//    G4cout
//       << "   Absorber: total energy: " << std::setw(7)
//                                        << G4BestUnit(fEnergyAbs,"Energy")
//       << "       total track length: " << std::setw(7)
//                                        << G4BestUnit(fTrackLAbs,"Length")
//       << G4endl
//       << "        Gap: total energy: " << std::setw(7)
//                                        << G4BestUnit(fEnergyGap,"Energy")
//       << "       total track length: " << std::setw(7)
//                                        << G4BestUnit(fTrackLGap,"Length")
//       << G4endl;
//  }
//}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
