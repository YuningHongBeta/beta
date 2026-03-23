// browse_calarr_cluster.C (PID from calarr version)  --- with cluster filtering ---
//
// 左：dE map（MeV） + (dE>=thr のbinに PIDラベル文字を重ね描き)
// 右：cluster ID map（任意）
// 上：文字サマリー（大きめ）
//
// 前提：rootファイルに TTree "calarr" があり、少なくとも
//   - eventID (int)
//   - dE_MeV  (std::vector<double>)  ※サイズ= nLayer*nSeg (通常225)
//   - pid     (std::vector<int>)     ※サイズ= nLayer*nSeg (通常225)
// が入っていること。
//
// 操作（標準入力）:
//   Enter       : 次のentry
//   b / back    : 前のentry
//   q / quit    : 終了
//   <number>    : entry番号へジャンプ（0..nent-1）
//   t <thr>     : threshold(MeV)変更（表示とクラスタ両方）
//   n <4|8>     : 近傍 4/8
//   w <0|1>     : wrapSeg ON/OFF（seg方向周期）
//   c <0|1>     : cluster map 表示 ON/OFF
//   z <val>     : dE map の Zmax 固定（val<=0 で解除）
//   log <0|1>   : 左パッド logz ON/OFF
//   txt <0|1>   : dE map TEXT ON/OFF（数字表示）
//   pid <0|1>   : PIDラベル表示 ON/OFF（文字表示）
//   tr <0|1>    : transpose ON/OFF
//   fy <0|1>    : flipY ON/OFF
//
// ★追加（今回）:
//   s <0|1>     : onlyOneCluster ON/OFF（ONなら nCluster==1 の entry だけ表示）
//   find1       : 次の nCluster==1 を探す（onlyOneClusterとは独立）
//   find0       : 次の nCluster==0 を探す
//
// 例:
//   root -l
//   root [0] .L browse_calarr_cluster.C
//   root [1] bcal("output/beta.root",15,15,0,0.2,8,true,true);

#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TPad.h"
#include "TPaveText.h"
#include "TLatex.h"

#include <iostream>
#include <string>
#include <vector>
#include <queue>
#include <algorithm>
#include <cstdio>   // ★ for snprintf
#include <cmath>

// --- helpers ---
struct RC { int r; int c; };

static const char* ZToSymbol(int Z)
{
  static const char* sym[] = {
    "",  "H","He",
    "Li","Be","B","C","N","O","F","Ne",
    "Na","Mg","Al","Si","P","S","Cl","Ar",
    "K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
    "Ga","Ge","As","Se","Br","Kr",
    "Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd",
    "In","Sn","Sb","Te","I","Xe",
    "Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu",
    "Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg",
    "Tl","Pb","Bi","Po","At","Rn",
    "Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr",
    "Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og"
  };
  if (Z < 1 || Z >= (int)(sizeof(sym)/sizeof(sym[0]))) return nullptr;
  return sym[Z];
}

static const char* PdgToShortLabel(int pdg)
{
  const int apdg = std::abs(pdg);

  // Ion PDG: 1000000000 + 10000*Z + 10*A + I
  if (apdg >= 1000000000) {
    const int Z = (apdg / 10000) % 1000;
    const int A = (apdg / 10) % 1000;
    const int I = apdg % 10;

    const char* sym = ZToSymbol(Z);
    static char buf[64];

    if (sym) {
      if (I == 0) std::snprintf(buf, sizeof(buf), "^{%d}%s", A, sym);
      else        std::snprintf(buf, sizeof(buf), "^{%d}%s*", A, sym);
    } else {
      std::snprintf(buf, sizeof(buf), "Z%dA%d", Z, A);
    }
    return buf;
  }

  switch (pdg) {
    case 22:   return "#gamma";
    case 11:   return "e-";
    case -11:  return "e+";
    case 13:   return "#mu-";
    case -13:  return "#mu+";
    case 211:  return "#pi+";
    case -211: return "#pi-";
    case 111:  return "#pi0";
    case 2212: return "p";
    case 2112: return "n";
    default:   return nullptr;
  }
}

static int CountHits_vec(const std::vector<double>& dE, int nRow, int nCol, double thrMeV)
{
  int nh = 0;
  const int nTot = nRow * nCol;
  for (int i = 0; i < std::min((int)dE.size(), nTot); i++) {
    if (dE[i] >= thrMeV) nh++;
  }
  return nh;
}

// connected components labeling (4 or 8 neighbors; seg wrap optional)
static int LabelClusters_vec(const std::vector<double>& dE,
                             int nRow, int nCol,
                             double thrMeV,
                             bool diag8,
                             bool wrapCol,
                             std::vector<std::vector<int>>& clusterID)
{
  clusterID.assign(nRow, std::vector<int>(nCol, 0));
  std::vector<std::vector<char>> vis(nRow, std::vector<char>(nCol, 0));

  auto isHit = [&](int r, int c) -> bool {
    const int idx = r * nCol + c;
    if (idx < 0 || idx >= (int)dE.size()) return false;
    return (dE[idx] >= thrMeV);
  };

  static const int dr4[4] = {-1, +1, 0, 0};
  static const int dc4[4] = {0, 0, -1, +1};
  static const int dr8[8] = {-1, -1, -1, 0, 0, +1, +1, +1};
  static const int dc8[8] = {-1, 0, +1, -1, +1, -1, 0, +1};

  int ncl = 0;

  for (int r = 0; r < nRow; r++) {
    for (int c = 0; c < nCol; c++) {
      if (!isHit(r, c) || vis[r][c]) continue;

      ncl++;
      std::queue<RC> q;
      vis[r][c] = 1;
      clusterID[r][c] = ncl;
      q.push({r, c});

      while (!q.empty()) {
        auto cur = q.front(); q.pop();

        const int nn = diag8 ? 8 : 4;
        for (int k = 0; k < nn; k++) {
          int nr = cur.r + (diag8 ? dr8[k] : dr4[k]);
          int nc = cur.c + (diag8 ? dc8[k] : dc4[k]);

          if (nr < 0 || nr >= nRow) continue;

          if (wrapCol) {
            if (nc < 0)      nc += nCol;
            if (nc >= nCol)  nc -= nCol;
          } else {
            if (nc < 0 || nc >= nCol) continue;
          }

          if (!vis[nr][nc] && isHit(nr, nc)) {
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

// Fill dE histogram (values below thr are forced to 0)
static void FillDE2D(TH2F* h,
                     const std::vector<double>& dE,
                     int nRow, int nCol,
                     double thrMeV,
                     bool transpose=false, bool flipY=false)
{
  h->Reset("ICES");

  for (int r = 0; r < nRow; r++) {
    for (int c = 0; c < nCol; c++) {
      const int idx = r * nCol + c;
      double v = (idx < (int)dE.size()) ? dE[idx] : 0.0;
      if (v < thrMeV) v = 0.0;

      int x = c + 1;
      int y = r + 1;
      if (flipY) y = nRow - r;
      if (transpose) std::swap(x, y);

      h->SetBinContent(x, y, v);
    }
  }
}

static void FillCluster2D(TH2F* hCl,
                          const std::vector<std::vector<int>>& clusterID,
                          int nRow, int nCol,
                          bool transpose=false, bool flipY=false)
{
  hCl->Reset("ICES");
  int maxID = 0;

  for (int r = 0; r < nRow; r++) {
    for (int c = 0; c < nCol; c++) {
      int id = clusterID[r][c];
      maxID = std::max(maxID, id);

      int x = c + 1;
      int y = r + 1;
      if (flipY) y = nRow - r;
      if (transpose) std::swap(x, y);

      hCl->SetBinContent(x, y, id);
    }
  }

  hCl->SetMinimum(0.0);
  hCl->SetMaximum(std::max(1, maxID));
}

// Draw PID labels on top of dE map (only bins where dE>=thr)
static void DrawPidOverlay(TH2F* hDE,
                           const std::vector<double>& dE,
                           const std::vector<int>& pid,
                           int nRow, int nCol,
                           double thrMeV,
                           bool transpose=false, bool flipY=false)
{
  TLatex lat;
  lat.SetTextAlign(22);
  lat.SetTextSize(0.030);
  lat.SetTextFont(42);
  lat.SetTextColor(kWhite);

  const int nTot = nRow * nCol;
  for (int r = 0; r < nRow; ++r) {
    for (int c = 0; c < nCol; ++c) {
      const int idx = r * nCol + c;
      if (idx >= (int)dE.size()) continue;
      if (idx >= (int)pid.size()) continue;
      if (idx >= nTot) continue;

      const double v = dE[idx];
      if (v < thrMeV) continue;

      const int pdg = pid[idx];
      if (pdg == 0) continue;

      int xbin = c + 1;
      int ybin = r + 1;
      if (flipY) ybin = nRow - r;
      if (transpose) std::swap(xbin, ybin);

      const double x = hDE->GetXaxis()->GetBinCenter(xbin);
      const double y = hDE->GetYaxis()->GetBinCenter(ybin);

      const char* lab = PdgToShortLabel(pdg);
      if (lab) lat.DrawLatex(x, y, lab);
      else     lat.DrawLatex(x, y, "N");
    }
  }
}

// -----------------------------------------------------------------------------
// bcal: browser
// -----------------------------------------------------------------------------
void bcal(const char* fname="output/beta.root",
          int nLayer=15, int nSeg=15,
          Long64_t start=0,
          double thrMeV=1.0,
          int neighbor=8,
          bool wrapSeg=true,
          bool showClusterMap=true,
          bool transpose=false, bool flipY=false,
          bool logz=false,
          bool drawText=false,
          double zMaxFixed=-1.0)
{
  gStyle->SetOptStat(0);

  if (!(neighbor == 4 || neighbor == 8)) {
    std::cout << "neighbor must be 4 or 8\n";
    return;
  }
  const int nTot = nLayer * nSeg;

  TFile* f = TFile::Open(fname, "READ");
  if (!f || f->IsZombie()) {
    std::cout << "Failed to open " << fname << "\n";
    return;
  }

  TTree* t = (TTree*)f->Get("calarr");
  if (!t) {
    std::cout << "No TTree named 'calarr'\n";
    return;
  }

  int eventID = -1;
  std::vector<double>* dE = nullptr;
  std::vector<int>*    pid = nullptr;

  t->SetBranchAddress("eventID", &eventID);
  t->SetBranchAddress("dE_MeV",  &dE);
  t->SetBranchAddress("pid",     &pid);

  const Long64_t nent = t->GetEntries();
  if (nent <= 0) {
    std::cout << "No entries\n";
    return;
  }

  Long64_t ev = start;
  if (ev < 0) ev = 0;
  if (ev >= nent) ev = nent - 1;

  // ---- UI state ----
  bool drawPID = true;

  // ★追加：cluster=1 のみ表示するか
  bool onlyOneCluster = false;

  // ---- Canvas / Pads ----
  TCanvas* c = new TCanvas("cCalArr", "calarr browser", 1300, 720);

  TPad* pSum = nullptr;
  TPad* pDE  = nullptr;
  TPad* pCL  = nullptr;
  TPaveText* sumBox = nullptr;

  auto rebuildPads = [&]() {
    c->Clear();

    pSum = new TPad("pSum", "pSum", 0.00, 0.88, 1.00, 1.00);
    pSum->SetFillColor(0);
    pSum->Draw();

    if (showClusterMap) {
      pDE = new TPad("pDE", "pDE", 0.00, 0.00, 0.60, 0.88);
      pCL = new TPad("pCL", "pCL", 0.60, 0.00, 1.00, 0.88);
      pDE->Draw(); pCL->Draw();
    } else {
      pDE = new TPad("pDE", "pDE", 0.00, 0.00, 1.00, 0.88);
      pCL = nullptr;
      pDE->Draw();
    }

    pSum->cd();
    if (sumBox) { delete sumBox; sumBox = nullptr; }
    sumBox = new TPaveText(0.01, 0.10, 0.99, 0.90, "NDC");
    sumBox->SetFillStyle(0);
    sumBox->SetBorderSize(1);
    sumBox->SetTextAlign(12);
    sumBox->SetTextSize(0.35);
    sumBox->Draw();
  };

  rebuildPads();

  // ---- Histograms ----
  auto getNx = [&]() { return transpose ? nLayer : nSeg; };
  auto getNy = [&]() { return transpose ? nSeg   : nLayer; };

  TH2F* hDE = new TH2F("hDE", "dE/PID map;segNo;layerNo",
                       getNx(), 0.0, (double)getNx(),
                       getNy(), 0.0, (double)getNy());
  hDE->SetMinimum(0.0);
  if (zMaxFixed > 0.0) hDE->SetMaximum(zMaxFixed);

  TH2F* hCl = nullptr;
  if (showClusterMap) {
    hCl = new TH2F("hCl", "cluster id;segNo;layerNo",
                   getNx(), 0.0, (double)getNx(),
                   getNy(), 0.0, (double)getNy());
  }

  std::vector<std::vector<int>> clusterID;

  // ---- helper: compute cluster count (and clusterID) for a given entry ----
  auto computeClustersForEntry = [&](Long64_t entry,
                                     int& outNh,
                                     int& outNcl,
                                     double& outSum,
                                     double& outMax) -> bool
  {
    if (entry < 0 || entry >= nent) return false;

    t->GetEntry(entry);
    if (!dE || !pid) return false;
    if ((int)dE->size() < nTot || (int)pid->size() < nTot) return false;

    const bool diag8 = (neighbor == 8);

    outSum = 0.0; outMax = 0.0;
    for (int i = 0; i < nTot; ++i) {
      outSum += (*dE)[i];
      outMax = std::max(outMax, (*dE)[i]);
    }

    outNh  = CountHits_vec(*dE, nLayer, nSeg, thrMeV);
    outNcl = LabelClusters_vec(*dE, nLayer, nSeg, thrMeV, diag8, wrapSeg, clusterID);
    return true;
  };

  // ---- helper: seek next entry with desired nCluster ----
  auto seekToNextMatching = [&](Long64_t from,
                                int dir,
                                int wantNcl,     // >=0 means exact match; <0 means "any"
                                bool useFilter) -> bool
  {
    Long64_t cur = from;
    while (true) {
      if (cur < 0 || cur >= nent) return false;

      int nh=0, ncl=0;
      double sum=0.0, mx=0.0;
      if (!computeClustersForEntry(cur, nh, ncl, sum, mx)) return false;

      if (!useFilter) { ev = cur; return true; }

      if (wantNcl >= 0) {
        if (ncl == wantNcl) { ev = cur; return true; }
      } else {
        // any
        ev = cur; return true;
      }

      cur += dir;
    }
  };

  // ---- redraw current ----
  auto redrawCurrent = [&]() {
    int nh=0, ncl=0;
    double sum=0.0, mx=0.0;

    if (!computeClustersForEntry(ev, nh, ncl, sum, mx)) {
      std::cout << "ERROR: failed to load entry " << ev << "\n";
      return;
    }

    // summary
    pSum->cd();
    sumBox->Clear();
    sumBox->AddText(Form("entry %lld / %lld   eventID=%d   sum(dE)=%.3f MeV   max=%.3f MeV",
                         ev, nent-1, eventID, sum, mx));
    sumBox->AddText(Form("thr=%.3f MeV   nhit=%d   nCluster=%d   neigh=%d   wrap=%d   only1=%d",
                         thrMeV, nh, ncl, neighbor, (int)wrapSeg, (int)onlyOneCluster));
    sumBox->AddText("keys: Enter(next) b(back) q(quit) <num>(jump)  t thr  n 4|8  w 0|1  c 0|1  z val  log 0|1  txt 0|1  pid 0|1  s 0|1  find1/find0");
    pSum->Modified(); pSum->Update();

    // dE map
    pDE->cd();
    pDE->SetLogz(logz ? 1 : 0);
    FillDE2D(hDE, *dE, nLayer, nSeg, thrMeV, transpose, flipY);
    hDE->Draw(drawText ? "COLZ TEXT" : "COLZ");

    if (drawPID) {
      DrawPidOverlay(hDE, *dE, *pid, nLayer, nSeg, thrMeV, transpose, flipY);
    }
    pDE->Modified(); pDE->Update();

    // cluster map
    if (showClusterMap && pCL && hCl) {
      pCL->cd();
      FillCluster2D(hCl, clusterID, nLayer, nSeg, transpose, flipY);
      hCl->Draw("COLZ");
      pCL->Modified(); pCL->Update();
    }

    c->cd();
    c->Modified(); c->Update();
    gSystem->ProcessEvents();
  };

  // ---- initial positioning if onlyOneCluster enabled? (starts off false) ----
  redrawCurrent();

  // interactive loop
  while (true) {
    std::cout << "[entry " << ev << "/" << (nent-1) << "] > ";
    std::string line;
    std::getline(std::cin, line);

    if (line == "q" || line == "quit") break;

    auto startsWith = [&](const std::string& s) {
      return line.size() >= s.size() && line.compare(0, s.size(), s) == 0;
    };

    // --- toggle onlyOneCluster ---
    if (startsWith("s ")) {
      try { onlyOneCluster = (std::stoi(line.substr(2)) != 0); } catch(...) {}
      if (onlyOneCluster) {
        // 現在位置から前後どちらでもいいが、とりあえず「現entryが合わないなら次の1へ」
        int nh=0, ncl=0; double sum=0, mx=0;
        computeClustersForEntry(ev, nh, ncl, sum, mx);
        if (ncl != 1) {
          if (!seekToNextMatching(ev, +1, 1, true)) {
            std::cout << "No further nCluster==1 entries.\n";
          }
        }
      }
      redrawCurrent();
      continue;
    }

    // --- find next exact cluster count (independent) ---
    if (line == "find1") {
      if (!seekToNextMatching(ev+1, +1, 1, true))
        std::cout << "No further nCluster==1 entries.\n";
      redrawCurrent();
      continue;
    }
    if (line == "find0") {
      if (!seekToNextMatching(ev+1, +1, 0, true))
        std::cout << "No further nCluster==0 entries.\n";
      redrawCurrent();
      continue;
    }

    // back
    if (line == "b" || line == "back") {
      if (ev > 0) ev--;
      if (onlyOneCluster) {
        if (!seekToNextMatching(ev, -1, 1, true)) {
          std::cout << "No previous nCluster==1 entries.\n";
        }
      }
      redrawCurrent();
      continue;
    }

    // Enter -> next
    if (line.size() == 0) {
      if (ev < nent - 1) ev++;
      else break;

      if (onlyOneCluster) {
        if (!seekToNextMatching(ev, +1, 1, true)) {
          std::cout << "No further nCluster==1 entries.\n";
          break;
        }
      }
      redrawCurrent();
      continue;
    }

    // thr
    if (startsWith("t ")) {
      try { thrMeV = std::stod(line.substr(2)); } catch(...) {}
      // filter mode なら “thr変更で cluster数も変わる” ので、現在から次の1へ合わせ直す
      if (onlyOneCluster) {
        int nh=0, ncl=0; double sum=0, mx=0;
        if (computeClustersForEntry(ev, nh, ncl, sum, mx) && ncl != 1) {
          if (!seekToNextMatching(ev, +1, 1, true))
            std::cout << "No nCluster==1 entries under this threshold.\n";
        }
      }
      redrawCurrent();
      continue;
    }

    // neighbor
    if (startsWith("n ")) {
      try {
        int nn = std::stoi(line.substr(2));
        if (nn == 4 || nn == 8) neighbor = nn;
      } catch(...) {}
      if (onlyOneCluster) {
        int nh=0, ncl=0; double sum=0, mx=0;
        if (computeClustersForEntry(ev, nh, ncl, sum, mx) && ncl != 1) {
          if (!seekToNextMatching(ev, +1, 1, true))
            std::cout << "No nCluster==1 entries under this neighbor setting.\n";
        }
      }
      redrawCurrent();
      continue;
    }

    // wrap
    if (startsWith("w ")) {
      try { wrapSeg = (std::stoi(line.substr(2)) != 0); } catch(...) {}
      if (onlyOneCluster) {
        int nh=0, ncl=0; double sum=0, mx=0;
        if (computeClustersForEntry(ev, nh, ncl, sum, mx) && ncl != 1) {
          if (!seekToNextMatching(ev, +1, 1, true))
            std::cout << "No nCluster==1 entries under this wrap setting.\n";
        }
      }
      redrawCurrent();
      continue;
    }

    // cluster map show/hide
    if (startsWith("c ")) {
      try { showClusterMap = (std::stoi(line.substr(2)) != 0); } catch(...) {}
      rebuildPads();
      if (showClusterMap && !hCl) {
        hCl = new TH2F("hCl", "cluster id;segNo;layerNo",
                       getNx(), 0.0, (double)getNx(),
                       getNy(), 0.0, (double)getNy());
      }
      redrawCurrent();
      continue;
    }

    // z max
    if (startsWith("z ")) {
      try {
        zMaxFixed = std::stod(line.substr(2));
        if (zMaxFixed > 0.0) hDE->SetMaximum(zMaxFixed);
      } catch(...) {}
      redrawCurrent();
      continue;
    }

    // log
    if (startsWith("log ")) {
      try { logz = (std::stoi(line.substr(4)) != 0); } catch(...) {}
      redrawCurrent();
      continue;
    }

    // txt
    if (startsWith("txt ")) {
      try { drawText = (std::stoi(line.substr(4)) != 0); } catch(...) {}
      redrawCurrent();
      continue;
    }

    // pid overlay
    if (startsWith("pid ")) {
      try { drawPID = (std::stoi(line.substr(4)) != 0); } catch(...) {}
      redrawCurrent();
      continue;
    }

    // transpose
    if (startsWith("tr ")) {
      try { transpose = (std::stoi(line.substr(3)) != 0); } catch(...) {}
      hDE->SetBins(getNx(), 0.0, (double)getNx(), getNy(), 0.0, (double)getNy());
      if (hCl) hCl->SetBins(getNx(), 0.0, (double)getNx(), getNy(), 0.0, (double)getNy());
      redrawCurrent();
      continue;
    }

    // flipY
    if (startsWith("fy ")) {
      try { flipY = (std::stoi(line.substr(3)) != 0); } catch(...) {}
      redrawCurrent();
      continue;
    }

    // number jump
    try {
      long long tgt = std::stoll(line);
      if (tgt < 0) tgt = 0;
      if (tgt >= nent) tgt = nent - 1;
      ev = (Long64_t)tgt;

      if (onlyOneCluster) {
        int nh=0, ncl=0; double sum=0, mx=0;
        if (computeClustersForEntry(ev, nh, ncl, sum, mx) && ncl != 1) {
          if (!seekToNextMatching(ev, +1, 1, true))
            std::cout << "No further nCluster==1 entries.\n";
        }
      }

      redrawCurrent();
    } catch (...) {
      // ignore
    }
  }

  std::cout << "Done.\n";
}