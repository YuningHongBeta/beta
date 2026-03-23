#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TStyle.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

void scan_target_peak(const char* fileName="output/beta.root",
                      double eMin=6.0,     // MeV
                      double eMax=10.0,     // MeV
                      int    maxPrint=50,  // 画面に出す最大件数
                      const char* outList="peak_events.txt")
{
  gStyle->SetOptStat(1110);

  TFile* f = TFile::Open(fileName, "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "ERROR: cannot open " << fileName << "\n";
    return;
  }

  // tree name: "target" 前提
  TTree* t = (TTree*)f->Get("target");
  if (!t) {
    std::cerr << "ERROR: TTree 'target' not found.\n";
    f->Close();
    return;
  }

  int eventID = -1;
  double dE_MeV = 0.0;

  // branch name: target:dE_MeV
  t->SetBranchAddress("eventID", &eventID);
  t->SetBranchAddress("dE_MeV",  &dE_MeV);

  auto* h = new TH1D("hTargetDE", "Target dE;Target dE (MeV);Events", 200, 0, 100);

  std::vector<int> peaks;
  const Long64_t nEnt = t->GetEntries();
  for (Long64_t i=0; i<nEnt; ++i) {
    t->GetEntry(i);
    h->Fill(dE_MeV);
    if (eMin <= dE_MeV && dE_MeV < eMax) peaks.push_back(eventID);
  }

  std::sort(peaks.begin(), peaks.end());
  peaks.erase(std::unique(peaks.begin(), peaks.end()), peaks.end());

  std::cout << "File: " << fileName << "\n";
  std::cout << "Entries(target): " << nEnt << "\n";
  std::cout << "Peak window: [" << eMin << ", " << eMax << "] MeV\n";
  std::cout << "N peak-like events: " << peaks.size() << "\n";

  // print some
  for (int i=0; i<(int)peaks.size() && i<maxPrint; ++i) {
    std::cout << "  eventID = " << peaks[i] << "\n";
  }
  if ((int)peaks.size() > maxPrint) std::cout << "  ...\n";

  // save list
  {
    std::ofstream ofs(outList);
    for (auto ev : peaks) ofs << ev << "\n";
  }
  std::cout << "Saved eventID list -> " << outList << "\n";

  // draw
  auto* c = new TCanvas("cTargetDE", "Target dE", 900, 600);
  h->Draw();
  c->Update();
}