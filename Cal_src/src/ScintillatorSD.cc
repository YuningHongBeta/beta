// ScintillatorSD.cc
#include "ScintillatorSD.hh"
#include "calEventAction.hh"
#include "calConstant.hh"

#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

#include "TROOT.h"
#include "TH1.h"
#include "TTree.h"
#include "TFile.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ScintillatorSD::ScintillatorSD(G4String name)
: G4VSensitiveDetector(name)
  //fHitsCollection(nullptr),fHCID(-1)
{
  collectionName.insert("ScintillatorColl");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ScintillatorSD::~ScintillatorSD()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ScintillatorSD::Initialize(G4HCofThisEvent* hce)
{
//  fHitsCollection 
//    = new ScintillatorHitsCollection(SensitiveDetectorName,collectionName[0]);
//  if (fHCID<0) {
//    fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection); 
//  }
//  hce->AddHitsCollection(fHCID,fHitsCollection);
//  
  etot_Scinti=0.;
  for (auto i=0;i<nb_cryst;i++) {
    //    fHitsCollection->insert(new ScintillatorHit(i));
    edepbuf_Scinti[i] =0.;
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool ScintillatorSD::ProcessHits(G4Step*step, G4TouchableHistory*)
{
  const G4StepPoint* preStepPoint = step->GetPreStepPoint();
  const G4StepPoint* postStepPoint = step->GetPostStepPoint();
  //G4String particleName = step->GetSecondaryInCurrentStep()->GetDefinition()->GetParticleName();
  //G4String proc = postStepPoint->GetProcessDefinedStep()->GetProcessName();
 
  auto edep = step->GetTotalEnergyDeposit();
  if (edep==0.) return true;
  
  auto touchable = step->GetPreStepPoint()->GetTouchable();
  auto detector = touchable->GetVolume();
  auto copyNo = detector->GetCopyNo();

  // get detector Name and ID 
  G4String detectorName = detector->GetName();    
  G4int            ID   = detector->GetCopyNo();  

  G4String particleName = step->GetTrack()->GetDefinition()->GetParticleName();


  edepbuf_Scinti[copyNo] += edep;  

  //auto hit = (*fHitsCollection)[copyNo];
  // check if it is first touch
  //if (!(hit->GetLogV())) {
    // fill volume information
    //hit->SetLogV(detector->GetLogicalVolume());
    //G4AffineTransform transform = touchable->GetHistory()->GetTopTransform();
    //transform.Invert();
    //hit->SetRot(transform.NetRotation());
    //hit->SetPos(transform.NetTranslation());

    //auto worldPos = preStepPoint->GetPosition();
    //auto GTime = preStepPoint->GetGlobalTime();    
    ////G4cout << copyNo << " " << GTime << " " << worldPos.x() << " " << worldPos.y() << " " << worldPos.z() << G4endl;
    //firsthitposXbuf[copyNo] = worldPos.x();  
    //firsthitposYbuf[copyNo] = worldPos.y();  
    //firsthitposZbuf[copyNo] = worldPos.z();  
    //firsthitTimebuf[copyNo] = GTime;  
    //getchar();
  //}


//  // add energy deposition
//  hit->AddEdep(edep);
  
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ScintillatorSD::EndOfEvent(G4HCofThisEvent* hce)
{
  //TH1D* Cal_dE= (TH1D*)gROOT-> FindObject("h_dE");
  for(int i=0;i<nb_cryst;i++){
    if(edepbuf_Scinti[i]>0){
      etot_Scinti+=edepbuf_Scinti[i];
      }
  }
  //G4cout << etot_Scinti << G4endl;
  //Cal_nhit->Fill(multi, 1);

  TTree* tree = (TTree*)gROOT-> FindObject("tree");
  tree->SetBranchAddress("dE_Scinti"	, edepbuf_Scinti);
  tree->SetBranchAddress("dEtot_Scinti"	, &etot_Scinti);

}

