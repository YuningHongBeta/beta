// browse_hitpat.C
#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TPad.h"
#include <iostream>
#include <string>
#include <cstdio>
#include <queue>
#include <vector>
#include <algorithm>
#include "TPaveText.h"
#include "TString.h"

struct RC
{
  int r;
  int c;
};

static int CountHits(const double *dE, int nRow, int nCol, double thrMeV)
{
  int nh = 0;
  for (int i = 0; i < nRow * nCol; i++)
    if (dE[i] >= thrMeV)
      nh++;
  return nh;
}

// clusterID[r][c] に 0(非hit) / 1..ncluster を入れる
static int LabelClusters(const double *dE,
                         int nRow, int nCol,
                         double thrMeV,
                         bool diag,
                         bool wrapCol,
                         std::vector<std::vector<int>> &clusterID)
{
  clusterID.assign(nRow, std::vector<int>(nCol, 0));

  std::vector<std::vector<char>> vis(nRow, std::vector<char>(nCol, 0));
  auto isHit = [&](int r, int c) -> bool
  {
    int idx = r * nCol + c;
    return (dE[idx] >= thrMeV);
  };

  static const int dr4[4] = {-1, +1, 0, 0};
  static const int dc4[4] = {0, 0, -1, +1};
  static const int dr8[8] = {-1, -1, -1, 0, 0, +1, +1, +1};
  static const int dc8[8] = {-1, 0, +1, -1, +1, -1, 0, +1};

  int ncl = 0;

  for (int r = 0; r < nRow; r++)
  {
    for (int c = 0; c < nCol; c++)
    {
      if (!isHit(r, c) || vis[r][c])
        continue;

      ncl++;
      std::queue<RC> q;
      vis[r][c] = 1;
      clusterID[r][c] = ncl;
      q.push({r, c});

      while (!q.empty())
      {
        auto cur = q.front();
        q.pop();

        const int nn = diag ? 8 : 4;
        for (int k = 0; k < nn; k++)
        {
          int nr = cur.r + (diag ? dr8[k] : dr4[k]);
          int nc = cur.c + (diag ? dc8[k] : dc4[k]);

          if (nr < 0 || nr >= nRow)
            continue;

          if (wrapCol)
          {
            if (nc < 0)
              nc += nCol;
            if (nc >= nCol)
              nc -= nCol;
          }
          else
          {
            if (nc < 0 || nc >= nCol)
              continue;
          }

          if (!vis[nr][nc] && isHit(nr, nc))
          {
            vis[nr][nc] = 1;
            clusterID[nr][nc] = ncl;
            q.push({nr, nc});
          }
        }
      }
    }
  }

  return ncl;
}

static void FillHit2D(TH2F *h, const double *dE,
                      int nRow, int nCol,
                      double thrMeV,
                      bool transpose = false, bool flipY = false)
{
  h->Reset("ICES");

  for (int j = 0; j < nRow; j++)
  {
    for (int ic = 0; ic < nCol; ic++)
    {
      int idx = j * nCol + ic; // idx = CopyNo想定
      double val = dE[idx];

      // threshold
      if (val < thrMeV)
        val = 0.0;

      int x = ic + 1;
      int y = j + 1;
      if (flipY)
        y = nRow - j;

      if (transpose)
      {
        int tx = y;
        int ty = x;
        x = tx;
        y = ty;
      }

      h->SetBinContent(x, y, val);
    }
  }
}

static void FillCluster2D(TH2F *hCl,
                          const std::vector<std::vector<int>> &clusterID,
                          int nRow, int nCol,
                          bool transpose = false, bool flipY = false)
{
  hCl->Reset("ICES");
  int maxID = 0;

  for (int j = 0; j < nRow; j++)
  {
    for (int ic = 0; ic < nCol; ic++)
    {
      int id = clusterID[j][ic];
      maxID = std::max(maxID, id);

      int x = ic + 1;
      int y = j + 1;
      if (flipY)
        y = nRow - j;

      if (transpose)
      {
        int tx = y;
        int ty = x;
        x = tx;
        y = ty;
      }

      hCl->SetBinContent(x, y, id);
    }
  }

  hCl->SetMinimum(0.0);
  hCl->SetMaximum(std::max(1, maxID));
}

void browse_hitpat(const char *fname = "test.root",
                   int nRow = 15, int nCol = 15,
                   Long64_t start = 0,
                   bool transpose = false, bool flipY = false,
                   bool logz = false,
                   // 追加：クラスタ定義
                   double thrMeV = 1.0,
                   bool diag = true,
                   bool wrapCol = true,
                   bool showClusterMap = true)
{
  gStyle->SetOptStat(0);

  TFile *f = TFile::Open(fname, "READ");
  if (!f || f->IsZombie())
  {
    printf("Failed to open %s\n", fname);
    return;
  }

  TTree *tree = (TTree *)f->Get("tree");
  if (!tree)
  {
    printf("No TTree named 'tree'\n");
    return;
  }

  const int nb = nRow * nCol;
  if (nb != 225)
  {
    printf("This macro currently expects 15x15 (=225). You gave %d.\n", nb);
    return;
  }

  double dE[225];
  tree->SetBranchAddress("dE", dE);

  Long64_t nent = tree->GetEntries();
  if (nent <= 0)
  {
    printf("No entries\n");
    return;
  }

  Long64_t ev = start;
  if (ev < 0)
    ev = 0;
  if (ev >= nent)
    ev = nent - 1;

  TCanvas *c = new TCanvas("cHit", "hit pattern browser", 1200, 650);

  // パッド構成：左=dE、右=clusterID（任意）
  TPad *p1 = nullptr;
  TPad *p2 = nullptr;

  if (showClusterMap)
  {
    p1 = new TPad("p1", "p1", 0.00, 0.00, 0.60, 1.00);
    p2 = new TPad("p2", "p2", 0.60, 0.00, 1.00, 1.00);
    p1->Draw();
    p2->Draw();
  }
  else
  {
    p1 = new TPad("p1", "p1", 0.00, 0.00, 1.00, 1.00);
    p1->Draw();
  }

  // --- dE map ---
  p1->cd();
  if (logz)
    p1->SetLogz();

  TH2F *h = new TH2F("hHit", ";col (icrys);row (layer)",
                     nCol, 0.5, nCol + 0.5,
                     nRow, 0.5, nRow + 0.5);
  h->SetMinimum(0.0);
  h->SetMaximum(100.0); // Z range fixed
  h->Draw("COLZ");

  // --- cluster map ---
  TH2F *hCl = nullptr;
  if (showClusterMap)
  {
    p2->cd();
    hCl = new TH2F("hCluster", "cluster id;col (icrys);row (layer)",
                   nCol, 0.5, nCol + 0.5,
                   nRow, 0.5, nRow + 0.5);
    hCl->Draw("COLZ");
  }

  std::vector<std::vector<int>> clusterID;

  while (true)
  {
    tree->GetEntry(ev);

    int nh = CountHits(dE, nRow, nCol, thrMeV);
    int ncl = LabelClusters(dE, nRow, nCol, thrMeV, diag, wrapCol, clusterID);

    // タイトル（dE map側にまとめて表示）
    h->SetTitle(Form("Event %lld | thr=%.1f MeV | nhit=%d | ncluster=%d  (Enter: next, b: back, q: quit, number: jump)",
                     ev, thrMeV, nh, ncl));

    p1->cd();
    FillHit2D(h, dE, nRow, nCol, thrMeV, transpose, flipY);
    p1->Modified();
    p1->Update();

    if (showClusterMap && hCl)
    {
      p2->cd();
      FillCluster2D(hCl, clusterID, nRow, nCol, transpose, flipY);
      p2->Modified();
      p2->Update();
    }

    c->cd();
    c->Modified();
    c->Update();
    gSystem->ProcessEvents();

    std::cout << "[event " << ev << "/" << (nent - 1) << "] > ";
    std::string line;
    std::getline(std::cin, line);

    if (line == "q" || line == "quit")
      break;

    if (line == "b" || line == "back")
    {
      if (ev > 0)
        ev--;
      continue;
    }

    if (line.size() == 0)
    {
      if (ev < nent - 1)
        ev++;
      else
        break;
      continue;
    }

    try
    {
      long long tgt = std::stoll(line);
      if (tgt < 0)
        tgt = 0;
      if (tgt >= nent)
        tgt = nent - 1;
      ev = (Long64_t)tgt;
    }
    catch (...)
    {
      // ignore unknown input
    }
  }

  std::cout << "Done.\n";
}

// --- 追加：ncluster を全イベントで数えて TH1I に詰める ---
void make_ncluster_hist(const char *fname = "test.root",
                        int nRow = 15, int nCol = 15,
                        double thrMeV = 1.0,
                        bool diag = false,
                        bool wrapCol = true,
                        int nbins = 20, int xmin = 0, int xmax = 20,
                        bool save = true,
                        const char *outname = "build/ncluster_hist.root")
{
  gStyle->SetOptStat(0);

  TFile *f = TFile::Open(fname, "READ");
  if (!f || f->IsZombie())
  {
    printf("Failed to open %s\n", fname);
    return;
  }

  TTree *tree = (TTree *)f->Get("tree");
  if (!tree)
  {
    printf("No TTree named 'tree'\n");
    return;
  }

  const int nb = nRow * nCol;
  if (nb != 225)
  {
    printf("This function currently expects 15x15 (=225). You gave %d.\n", nb);
    return;
  }

  double dE[225];
  tree->SetBranchAddress("dE", dE);

  Long64_t nent = tree->GetEntries();
  if (nent <= 0)
  {
    printf("No entries\n");
    return;
  }

  TH1F *hNcl = new TH1F("hNcluster",
                        Form("ncluster (thr=%.1f MeV, %s, wrapCol=%d);ncluster;events",
                             thrMeV, diag ? "8-neigh" : "4-neigh", (int)wrapCol),
                        nbins, xmin, xmax);

  std::vector<std::vector<int>> clusterID;

  for (Long64_t ev = 0; ev < nent; ev++)
  {
    tree->GetEntry(ev);
    int ncl = LabelClusters(dE, nRow, nCol, thrMeV, diag, wrapCol, clusterID);
    hNcl->Fill(ncl);
  }

  TCanvas *c = new TCanvas("cNcluster", "ncluster histogram", 750, 900);
  hNcl->Draw();

  // ---- ここから統計をCanvas上に表示 ----
  Int_t nTotal = (Int_t)hNcl->GetEntries();

  // 注意：このヒストは xmin-0.5 なので bin1が ncluster=xmin に対応
  // xmin=0 のとき bin1=0 cluster, bin2=1 cluster
  Int_t binZero = hNcl->FindBin(0);
  Int_t binSingle = hNcl->FindBin(1);

  Int_t nZero = (Int_t)hNcl->GetBinContent(binZero);
  Int_t nSingle = (Int_t)hNcl->GetBinContent(binSingle);
  Int_t nMulti = nTotal - nZero - nSingle;

  double pZero = 0, pSingle = 0, pMulti = 0;
  if (nTotal > 0)
  {
    pZero = 100.0 * nZero / nTotal;
    pSingle = 100.0 * nSingle / nTotal;
    pMulti = 100.0 * nMulti / nTotal;
  }

  // 右上にテキストボックス
  TPaveText *pt = new TPaveText(0.50, 0.50, 0.88, 0.88, "NDC");
  pt->SetFillStyle(0);
  pt->SetBorderSize(1);
  pt->SetTextAlign(12); // left, center vertically
  pt->SetTextSize(0.025);

  pt->AddText(Form("Total events: %d", nTotal));
  pt->AddText(Form("Zero-cluster:  %d (%.2f%%)", nZero, pZero));
  pt->AddText(Form("Single-cluster:%d (%.2f%%)", nSingle, pSingle));
  pt->AddText(Form("Multi-cluster: %d (%.2f%%)", nMulti, pMulti));

  // pt->Draw();
  c->Modified();
  c->Update();

  if (save)
  {
    TFile *fo = TFile::Open(outname, "RECREATE");
    if (fo && !fo->IsZombie())
    {
      hNcl->Write();
      fo->Close();
    }
  }
}