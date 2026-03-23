void drawMean(char *fname){
  int nb_cryst=50; // number of Layers
  double absoThickness = 50;//[cm]
  double LayerThickness = absoThickness/double(nb_cryst);
  double depth_dE;
  double depth_Hit;

  int eveno;
  int multiHit;
  int cellNoFlag[200];
  double dEtot;
  double dE[200];
  double dE_alleve[200];
  double dE_single;

  TFile *ifile = TFile::Open(fname);
  TTree *tree = (TTree*)ifile->Get("tree");
  tree->SetBranchAddress("eveno"      ,&eveno);
  tree->SetBranchAddress("multiHit"   ,&multiHit);
  tree->SetBranchAddress("cellNoFlag"   ,cellNoFlag);
  tree->SetBranchAddress("dEtot"      ,&dEtot);
  tree->SetBranchAddress("dE"      ,dE);


  double xmin = 0;
  double xmax = absoThickness;
  int xbin = (xmax -xmin)*1;
  TH1F *h_multi = new TH1F("h_multi","h_multi",xbin ,xmin, xmax);
  TGraph* g = new TGraph();

//  int event_max = tree->GetEntries();
//  for(int ev=0; ev<event_max; ev++){
//    dE_single = 0.;
//    tree->GetEntry(ev);
//
//    for(int id=0; id<nb_cryst; id++){
//      dE_alleve[id] = dE_alleve[id] + dE[id];
//    }
//
//    depth_Hit = 0;
//    depth_Hit = LayerThickness*multiHit + LayerThickness/2.;
//    h_multi->Fill(depth_Hit);
//  } //for ev                                                   

  TH1F *h[nb_cryst];
  for(int id=0; id<nb_cryst; id++){
    h[id] = new TH1F(Form("h%d",id),Form("h%d",id), 100, 0, 100);
  }

  double sum;
  double mean;
  for(int id=0; id<nb_cryst; id++){
    tree->Project(Form("h%d",id),Form("dE[%d]",id));
    mean = h[id]->GetMean();
    //cout << "##########" << mean << endl;
    //depth_dE = LayerThickness*id + LayerThickness/2.;
    g->SetPoint(id,id,mean);
    //g->SetPoint(id,depth_dE,dE_alleve[id]);
  }//for id

#if 1
  TCanvas *c = new TCanvas("c","c",700,400);
  gStyle->SetOptStat(0);
  //c->Divide(2,1);

  c->cd(1);
  gStyle->SetOptStat(0);
  g->SetTitle("dE (mean)");
  g->GetXaxis()->SetTitle("depth [cm]");
  g->GetYaxis()->SetTitle("dE_mean");

  g->Draw("apl");

//  c->cd(2);
//  gStyle->SetOptStat(0);
//  h_multi->SetTitle("multi Hit");
//  h_multi->GetXaxis()->SetTitle("depth [cm]");
//  h_multi->GetYaxis()->SetTitle("counts");
//  h_multi->Draw();

#endif

}
void drawBoxdE(char *fname){
  int nb_cryst=50; // number of Layers
  double absoThickness = 50;//[cm]
  double LayerThickness = absoThickness/double(nb_cryst);
  double depth_dE;
  double depth_Hit;

  int eveno;
  int multiHit;
  int cellNoFlag[200];
  double dEtot;
  double dE[200];
  double dE_alleve[200];
  double dE_single;

  TFile *ifile = TFile::Open(fname);
  TTree *tree = (TTree*)ifile->Get("tree");
  tree->SetBranchAddress("eveno"      ,&eveno);
  tree->SetBranchAddress("multiHit"   ,&multiHit);
  tree->SetBranchAddress("cellNoFlag"   ,cellNoFlag);
  tree->SetBranchAddress("dEtot"      ,&dEtot);
  tree->SetBranchAddress("dE"      ,dE);


  double xmin = 0;
  double xmax = absoThickness;
  int xbin = (xmax -xmin)*1;
  TH1F *h_multi = new TH1F("h_multi","h_multi",xbin ,xmin, xmax);
  TGraph* g = new TGraph();

  int event_max = tree->GetEntries();
  for(int ev=0; ev<event_max; ev++){
    dE_single = 0.;
    tree->GetEntry(ev);

    for(int id=0; id<nb_cryst; id++){
      dE_alleve[id] = dE_alleve[id] + dE[id];
    }

    depth_Hit = 0;
    depth_Hit = LayerThickness*multiHit + LayerThickness/2.;
    h_multi->Fill(depth_Hit);
  } //for ev                                                                                                                          
  depth_dE = 0;
  double sum;
  for(int id=0; id<nb_cryst; id++){
    depth_dE = LayerThickness*id + LayerThickness/2.;
    dE_alleve[id] = dE_alleve[id]/event_max;
    sum = sum + dE_alleve[id];
    g->SetPoint(id,depth_dE,dE_alleve[id]);
  }//for id


  TCanvas *c = new TCanvas("c","c",700,400);
  gStyle->SetOptStat(0);
  c->Divide(2,1);

  c->cd(1);
  gStyle->SetOptStat(0);
  g->SetTitle("dE (average)");
  g->GetXaxis()->SetTitle("depth [cm]");
  g->GetYaxis()->SetTitle("dE_average");

  g->Draw("apl");

  c->cd(2);
  gStyle->SetOptStat(0);
  h_multi->SetTitle("multi Hit");
  h_multi->GetXaxis()->SetTitle("depth [cm]");
  h_multi->GetYaxis()->SetTitle("counts");
  h_multi->Draw();



}
void drawdE(char *fname){
  const int nLayer=13;
  const int sector=8;
  const int nb_cryst=nLayer*sector;
  int eveno;
  int multiHit;
  int cellNoFlag[200];
  double dEtot;
  double dE[200];
  double dE_single;

  TFile *ifile = TFile::Open(fname);
  TTree *tree = (TTree*)ifile->Get("tree");
  tree->SetBranchAddress("eveno"      ,&eveno);
  tree->SetBranchAddress("multiHit"   ,&multiHit);
  tree->SetBranchAddress("cellNoFlag"   ,cellNoFlag);
  tree->SetBranchAddress("dEtot"      ,&dEtot);
  tree->SetBranchAddress("dE"      ,dE);


  int xmin = 0;
  int xmax = 200;
  int xbin = (xmax -xmin)*1;
  int ymin = 0;
  int ymax = 200;
  int ybin = (ymax - ymin)/1;
  TH2F *h = new TH2F("h","h",xbin ,xmin, xmax, ybin, ymin, ymax);
  TH2F *h2 = new TH2F("h2","h2",xbin ,xmin, xmax, ybin, ymin, ymax);

  TH1F *h_dE = new TH1F("h_dE","h_dE",xbin ,xmin, xmax);

  //int event_max = 1; 
  int event_max = tree->GetEntries();
  
  int layer = 0;
  int cell=1;


  for(int ev=0; ev<event_max; ev++){
    dE_single = 0.;
    tree->GetEntry(ev);

    for(int id=0; id<nb_cryst; id++){
      layer = id/sector;

      if(dE[id]>0){
	cout << layer << " " << cell << " " << dE[id] << endl;
	if(dE[id]>dE_single) { 
	  dE_single = dE[id];
	}

	cell++;
	if(cell==sector+1){  cell =1;    }
      }//if(dE[id]>0)

    }
    cout << "KORE!!" << dE_single << endl;
    h_dE->Fill(dE_single);
    getchar();


  } //for ev                                                                                                                          

  TCanvas *c = new TCanvas("c","c",400,400);
  gStyle->SetOptStat(0);
  h_dE->Draw();


}
