// CalorimeterSD.cc
#include "CalorimeterSD.hh"
#include "CalorimeterHit.hh"
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

CalorimeterSD::CalorimeterSD(G4String name)
: G4VSensitiveDetector(name), 
  fHitsCollection(nullptr), fHCID(-1)
{
  collectionName.insert("calorimeterColl");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CalorimeterSD::~CalorimeterSD()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void CalorimeterSD::Initialize(G4HCofThisEvent* hce)
{
  fHitsCollection 
    = new CalorimeterHitsCollection(SensitiveDetectorName,collectionName[0]);
  if (fHCID<0) {
    fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection); 
  }
  hce->AddHitsCollection(fHCID,fHitsCollection);
  
  for (auto i=0;i<nb_cryst;i++) {
    fHitsCollection->insert(new CalorimeterHit(i));
  }

  multiHit=0;  
  etot=0.;  
  for (auto i=0;i<nb_cryst;i++) {
    edepbuf[i] =0.;
    cellNoFlag[i] = 0;
    firsthitposXbuf[i] = 0.;
    firsthitposYbuf[i] = 0.;
    firsthitposZbuf[i] = 0.;
    firsthitTimebuf[i] = 0.;
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool CalorimeterSD::ProcessHits(G4Step*step, G4TouchableHistory*)
{
  const G4StepPoint* preStepPoint = step->GetPreStepPoint();
  const G4StepPoint* postStepPoint = step->GetPostStepPoint();
  //G4String particleName = step->GetSecondaryInCurrentStep()->GetDefinition()->GetParticleName();
  //G4String proc = postStepPoint->GetProcessDefinedStep()->GetProcessName();
 
  //  const std::vector<const G4Track*>* secondary = step->GetSecondaryInCurrentStep();
  //  for(int lp=0; lp<(*secondary).size();lp++){
  //    auto particle = (*secondary)[lp]->GetDefinition();
  //    G4String pName = (*secondary)[lp]->GetDefinition()->GetParticleName();
  //    auto touchable = step->GetPreStepPoint()->GetTouchable();
  //    auto detector = touchable->GetVolume();
  //    auto copyNo = detector->GetCopyNo();
  //
  //  }
  

  auto edep = step->GetTotalEnergyDeposit();
  if (edep==0.) return true;
  
  auto touchable = step->GetPreStepPoint()->GetTouchable();
  auto detector = touchable->GetVolume();
  auto copyNo = detector->GetCopyNo();

  // get detector Name and ID 
  G4String detectorName = detector->GetName();    
  G4int            ID   = detector->GetCopyNo();  

  G4String particleName = step->GetTrack()->GetDefinition()->GetParticleName();


  edepbuf[copyNo] += edep;  

  auto hit = (*fHitsCollection)[copyNo];
  // check if it is first touch
  if (!(hit->GetLogV())) {
    // fill volume information
    hit->SetLogV(detector->GetLogicalVolume());
    //G4AffineTransform transform = touchable->GetHistory()->GetTopTransform();
    //transform.Invert();
    //hit->SetRot(transform.NetRotation());
    //hit->SetPos(transform.NetTranslation());

    auto worldPos = preStepPoint->GetPosition();
    auto GTime = preStepPoint->GetGlobalTime();    
    //G4cout << copyNo << " " << GTime << " " << worldPos.x() << " " << worldPos.y() << " " << worldPos.z() << G4endl;
    firsthitposXbuf[copyNo] = worldPos.x();  
    firsthitposYbuf[copyNo] = worldPos.y();  
    firsthitposZbuf[copyNo] = worldPos.z();  
    firsthitTimebuf[copyNo] = GTime;  
    //getchar();
  }


//  // add energy deposition
//  hit->AddEdep(edep);
  
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void CalorimeterSD::EndOfEvent(G4HCofThisEvent* hce)
{

  ///TH1D* Cal_hit= (TH1D*)gROOT-> FindObject("g1");
  ///TH1D* Cal_nhit= (TH1D*)gROOT-> FindObject("g2");
  TH1D* Cal_dE= (TH1D*)gROOT-> FindObject("h_dE");

  multiHit=0;

  for(int i=0;i<nb_cryst;i++){

    cellNoFlag[i] = 0;
    if(edepbuf[i]>0){
      //Cal_hit->Fill(i, 1);
      //Cal_tedep->Fill(i,edepbuf[i]);
      cellNoFlag[i] = 1;
      etot+=edepbuf[i];
      multiHit++;
      }
  }

  //Cal_nhit->Fill(multi, 1);
  Cal_dE->Fill(etot);

  TTree* tree = (TTree*)gROOT-> FindObject("tree");
  //tree->SetBranchAddress("cellNoFlag"	, cellNoFlag);
  tree->SetBranchAddress("dE"		, edepbuf);
  tree->SetBranchAddress("dEtot_Cal"	, &etot);
  tree->SetBranchAddress("multiHit"	, &multiHit);
  //tree->SetBranchAddress("FirstHitPosX"		, firsthitposXbuf);
  //tree->SetBranchAddress("FirstHitPosY"		, firsthitposYbuf);
  //tree->SetBranchAddress("FirstHitPosZ"		, firsthitposZbuf);
  //tree->SetBranchAddress("FirstHitTime"		, firsthitTimebuf);

}

