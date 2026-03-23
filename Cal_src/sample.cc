const int nLayer = 15; // 30
const int sector = 15; // 60
const int nb_cryst = nLayer * sector;
int eveno;
double FirstHitTime[nb_cryst];
double FirstHitPosX[nb_cryst];
double FirstHitPosY[nb_cryst];
double FirstHitPosZ[nb_cryst];
double dE_Cal[nb_cryst];
double dE_Scinti[nb_cryst];
double dEtot_Cal;

void HitPat()
{
  // Draw hit pattern

  // open root file and get tree, branch
  TFile *ifile = TFile::Open("../build/test.root");
  // TFile *ifile = TFile::Open(fname);
  TTree *tree = (TTree *)ifile->Get("tree");
  tree->SetBranchAddress("eveno", &eveno);
  tree->SetBranchAddress("dE", dE_Cal);
  tree->SetBranchAddress("dE_Scinti", dE_Scinti);
  tree->SetBranchAddress("dEtot_Cal", &dEtot_Cal);

  TCanvas *c = new TCanvas("c", "c");
  gStyle->SetOptStat(0);

  // difine histogram
  double xmin = 10;
  double xmax = 200;
  int xbin = (xmax - xmin) * 1;
  TH1F *h[nLayer][sector];
  for (int k = 0; k < nLayer; k++)
  {
    for (int l = 0; l < sector; l++)
    {
      h[k][l] = new TH1F(Form("h_%d_%d", k, l), Form("h_%d_%d", k, l), xbin, xmin, xmax);
    }
  }

  TH2F *h_hitPat = new TH2F("h_hitPat", "h_hitPat", nLayer, 0, nLayer, sector, 0, sector);
  // auto frame = c->DrawFrame(0, 0, nLayer, sector);

#if 1
  int event_max = tree->GetEntries();
  event_max = 1000;

  for (int i = 0; i < event_max; i++)
  {
    // h_hitPat->Reset();
    tree->GetEntry(i);
    if (dEtot_Cal > 100)
    {
      continue;
    }

    int lay = 0;
    int sec = 0;
    for (int id = 0; id < nb_cryst; id++)
    {
      h[lay][sec]->Fill(dE_Cal[id]);
      h_hitPat->Fill(lay, sec, dE_Cal[id]);
      if (dE_Cal[id] > 0)
      {
        // cout  << lay << " " << sec << " " << dE[id] <<  "  " << id << endl;
      }
      sec++;
      if (sec == sector)
      {
        lay++;
        sec = 0;
      }

    } // for id
    cout << i << "  dEtot : " << dEtot_Cal << endl;

    gPad->SetLogz();
    h_hitPat->Draw("colz");
    c->Update();
    // getchar();
    // break;

  } // for i
#endif

#if 0
  //draw   ex. 0Layer 0sectorのdE分布を描く
  TCanvas *c1 = new TCanvas("c1","c1");
  gStyle->SetOptStat(0);
  c1->cd(1);
  gStyle->SetOptStat(0);
  h[0][0]->SetTitle("");
  h[0][0]->GetXaxis()->SetTitle("dE [MeV]");
  h[0][0]->GetYaxis()->SetTitle("counts");
  h[0][0]->Draw();
#endif
}
void angle(int num)
{
  TVector3 gamma1;
  TVector3 gamma2;

  // open root file and get tree, branch
  TFile *ifile = TFile::Open("../Cal-build/test.root");
  // TFile *ifile = TFile::Open(fname);
  TTree *tree = (TTree *)ifile->Get("tree");
  tree->SetBranchAddress("eveno", &eveno);
  tree->SetBranchAddress("dE", dE_Cal);
  tree->SetBranchAddress("dE_Scinti", dE_Scinti);
  tree->SetBranchAddress("dEtot_Cal", &dEtot_Cal);
  tree->SetBranchAddress("FirstHitTime", FirstHitTime);
  //  tree->SetBranchAddress("FirstHitPosX"      ,FirstHitPosX);
  //  tree->SetBranchAddress("FirstHitPosY"      ,FirstHitPosY);
  //  tree->SetBranchAddress("FirstHitPosZ"      ,FirstHitPosZ);

  TCanvas *c = new TCanvas("c", "c");
  gStyle->SetOptStat(0);
  c->Divide(1, 2);

  // difine histogram
  TH2F *h_hitPat = new TH2F("h_hitPat", "h_hitPat", nLayer, 0, nLayer, sector, 0, sector);
  TH2F *h_hitPat2 = new TH2F("h_hitPat2", "h_hitPat2", nLayer, 0, nLayer, sector, 0, sector);

  tree->GetEntry(num);
  // if(dEtot > 100){ continue; }

  int lay = 0;
  int sec = 0;
  for (int id = 0; id < nb_cryst; id++)
  {
    h_hitPat->Fill(lay, sec, dE_Cal[id]);
    if (dE_Cal[id] > 0)
    {
      // if(dE[id]>0 && FirstHitTime[id] < 2.0){
      cout << id << " " << lay << " " << sec << " " << dE_Cal[id] << " "
           << FirstHitTime[id] << " "
           //<< FirstHitPosX[id] << " " << FirstHitPosY[id] << " " << FirstHitPosZ[id] << " "
           << endl;
      h_hitPat2->Fill(lay, sec, dE_Scinti[id]);
    }
    sec++;
    if (sec == sector)
    {
      lay++;
      sec = 0;
    }

  } // for id
  cout << num << "  dEtot : " << dEtot_Cal << endl;

#if 0
  int id1 = 12;
  int id2 = 764;
  gamma1.SetXYZ(FirstHitPosX[id1], FirstHitPosY[id1], FirstHitPosZ[id1]);
  gamma2.SetXYZ(FirstHitPosX[id2], FirstHitPosY[id2], FirstHitPosZ[id2]);
  Double_t angle = gamma1.Angle(gamma2);
  cout << "angle : " << angle*180/3.14 << "[deg]" << endl;
  //getchar();
#endif

  c->cd(1);
  gPad->SetLogz();
  h_hitPat->Draw("colz");

  c->cd(2);
  gPad->SetLogz();
  h_hitPat2->Draw("colz");

  c->Update();
}
