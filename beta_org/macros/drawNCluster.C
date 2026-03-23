#include <TFile.h>
#include <TTree.h>
#include <TH1I.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TPaveText.h>
#include <TH1D.h>
#include <TLegend.h>

#include <vector>
#include <queue>
#include <iostream>
#include <unordered_map>
#include <cmath>
#include <algorithm> // sort

enum ClCat
{
  kClCat_Neutron = 1,
  kClCat_Gamma = 2,
  kClCat_Elec = 3,
  kClCat_Proton = 4,
  kClCat_Muon = 5,
  kClCat_Pion = 6,
  kClCat_Ion = 7,
  kClCat_Other = 8
};

static inline int idx2d(int layer, int seg, int nSeg)
{
  return layer * nSeg + seg;
}

static inline std::uint64_t keyEvCopy(int ev, int copyNo)
{
  return ((std::uint64_t)(std::uint32_t)ev << 32) | (std::uint32_t)copyNo;
}

static int pdgCategory(int pdg)
{
  if (pdg == 2112)
    return kClCat_Neutron;
  if (pdg == 22)
    return kClCat_Gamma;
  if (pdg == 11 || pdg == -11)
    return kClCat_Elec;
  if (pdg == 2212)
    return kClCat_Proton;
  if (pdg == 13 || pdg == -13)
    return kClCat_Muon;
  if (pdg == 211 || pdg == -211)
    return kClCat_Pion;
  if (std::abs(pdg) > 1000000000)
    return kClCat_Ion;
  return kClCat_Other;
}

static const char *catName(int c)
{
  switch (c)
  {
  case kClCat_Neutron:
    return "n";
  case kClCat_Gamma:
    return "gamma";
  case kClCat_Elec:
    return "e±";
  case kClCat_Proton:
    return "p";
  case kClCat_Muon:
    return "mu±";
  case kClCat_Pion:
    return "pi±";
  case kClCat_Ion:
    return "ion";
  default:
    return "other";
  }
}

struct ClInfo
{
  int size = 0;
  double sumE = 0.0;
};

static int countClustersWithThr(const std::vector<double> &de,
                                int nLayer, int nSeg,
                                double thrHitMeV,
                                double thrClMeV, // <0 なら無効
                                bool wrapSeg,
                                int neighbor,
                                int *firstClSize = nullptr,
                                double *firstClSumE = nullptr)
{
  const int nTot = nLayer * nSeg;

  std::vector<char> hit(nTot, 0), vis(nTot, 0);
  for (int k = 0; k < nTot; ++k)
  {
    if (de[k] > thrHitMeV)
      hit[k] = 1;
  }

  int nCl = 0;
  int firstSize = 0;
  double firstSum = 0.0;

  for (int layer = 0; layer < nLayer; ++layer)
  {
    for (int seg = 0; seg < nSeg; ++seg)
    {
      const int start = idx2d(layer, seg, nSeg);
      if (!hit[start] || vis[start])
        continue;

      // 1クラスタ分 BFS
      int thisSize = 0;
      double thisSum = 0.0;

      std::queue<std::pair<int, int>> q;
      q.push({layer, seg});
      vis[start] = 1;

      while (!q.empty())
      {
        auto [ly, sg] = q.front();
        q.pop();

        const int idx = idx2d(ly, sg, nSeg);
        ++thisSize;
        thisSum += de[idx];

        for (int dly = -1; dly <= 1; ++dly)
        {
          const int nly = ly + dly;
          if (nly < 0 || nly >= nLayer)
            continue;

          for (int dsg = -1; dsg <= 1; ++dsg)
          {
            if (dly == 0 && dsg == 0)
              continue;
            if (neighbor == 4 && (dly != 0 && dsg != 0))
              continue;

            int nsg = sg + dsg;
            if (wrapSeg)
            {
              if (nsg >= nSeg)
                nsg -= nSeg;
              if (nsg < 0)
                nsg += nSeg;
            }
            else
            {
              if (nsg < 0 || nsg >= nSeg)
                continue;
            }

            const int nb = idx2d(nly, nsg, nSeg);
            if (hit[nb] && !vis[nb])
            {
              vis[nb] = 1;
              q.push({nly, nsg});
            }
          }
        }
      }

      // cluster threshold（総和で判定）
      if (thrClMeV >= 0.0 && thisSum <= thrClMeV)
      {
        // このクラスタは「無効」なので数えない
        continue;
      }

      // 採用
      ++nCl;
      if (nCl == 1)
      {
        firstSize = thisSize;
        firstSum = thisSum;
      }
    }
  }

  if (firstClSize)
    *firstClSize = firstSize;
  if (firstClSumE)
    *firstClSumE = firstSum;
  return nCl;
}

void drawNClusterSummary(const char *fileName = "output/beta.root",
                         int nLayer = 15,
                         int nSeg = 15,
                         double thrMeV = 1.0,
                         double thrClMeV = 1.0,
                         bool wrapSeg = true,
                         int neighbor = 4,   // <-- 4 or 8
                         bool overlay = true,
                         double yMax = 0.0) // <=0 means auto
{
  if (!overlay && neighbor != 4 && neighbor != 8)
  {
    std::cerr << "ERROR: neighbor must be 4 or 8 (given " << neighbor << ")\n";
    return;
  }

  const int nTot = nLayer * nSeg;

  TFile *f = TFile::Open(fileName, "READ");
  if (!f || f->IsZombie())
  {
    std::cerr << "ERROR: cannot open file: " << fileName << "\n";
    return;
  }

  TTree *t = (TTree *)f->Get("calarr");
  if (!t)
  {
    std::cerr << "ERROR: TTree 'calarr' not found.\n";
    f->Close();
    return;
  }

  int ev = -1;
  std::vector<double> *de = nullptr;
  t->SetBranchAddress("eventID", &ev);
  t->SetBranchAddress("dE_MeV", &de);

  struct SummaryStats
  {
    Long64_t nEvents = 0;
    Long64_t n0 = 0;
    Long64_t n1 = 0;
    Long64_t nMulti = 0;
  };

  auto fillHist = [&](int neigh, TH1I *h, SummaryStats &st)
  {
    const Long64_t nEnt = t->GetEntries();
    for (Long64_t i = 0; i < nEnt; ++i)
    {
      t->GetEntry(i);
      if (!de)
        continue;
      if ((int)de->size() < nTot)
        continue;

      int nCluster = countClustersWithThr(*de, nLayer, nSeg,
                                          thrMeV, thrClMeV,
                                          wrapSeg, neigh);

      h->Fill(nCluster);
      ++st.nEvents;
      if (nCluster == 0)
        ++st.n0;
      else if (nCluster == 1)
        ++st.n1;
      else
        ++st.nMulti;
    }
  };

  gStyle->SetOptStat(0);

  auto *c = new TCanvas("cNClusterSummary",
                        "nCluster summary",
                        750, 900);
  c->SetLeftMargin(0.20);
  c->SetRightMargin(0.06);
  c->SetTopMargin(0.08);
  c->SetBottomMargin(0.12);

  if (!overlay)
  {
    const char *neighLabel = (neighbor == 4) ? "4-neighbor" : "8-neighbor";
    auto *h = new TH1I(Form("hNCluster_%d", neighbor),
                       Form("nCluster (%s, thr=%.3f MeV);nCluster;Events", neighLabel, thrMeV),
                       20, 0, 20);
    SummaryStats st;
    fillHist(neighbor, h, st);

    if (st.nEvents == 0)
    {
      std::cerr << "ERROR: no valid events found.\n";
      f->Close();
      return;
    }

    h->SetLineWidth(2);
    h->GetYaxis()->SetTitleOffset(2.0);
    h->GetXaxis()->SetTitleOffset(1.1);
    if (yMax > 0.0)
      h->SetMaximum(yMax);
    h->Draw();

    const double r0 = 100.0 * (double)st.n0 / (double)st.nEvents;
    const double r1 = 100.0 * (double)st.n1 / (double)st.nEvents;
    const double rM = 100.0 * (double)st.nMulti / (double)st.nEvents;

    auto *p = new TPaveText(0.55, 0.65, 0.88, 0.88, "NDC");
    p->SetFillColor(0);
    p->SetTextAlign(12);
    p->AddText(Form("thr: %.3f MeV", thrMeV));
    p->AddText(Form("cluster thr: %.3f MeV", thrClMeV));
    p->AddText(Form("Events: %lld", st.nEvents));
    p->AddText(Form("no cluster (0): %lld  (%.2f%%)", st.n0, r0));
    p->AddText(Form("one cluster (1): %lld  (%.2f%%)", st.n1, r1));
    p->AddText(Form("multi (>=2): %lld  (%.2f%%)", st.nMulti, rM));
    p->Draw();
  }
  else
  {
    auto *h4 = new TH1I("hNCluster_4",
                        Form("nCluster (thr=%.3f MeV);nCluster;Events", thrMeV),
                        20, 0, 20);
    auto *h8 = new TH1I("hNCluster_8",
                        Form("nCluster (thr=%.3f MeV);nCluster;Events", thrMeV),
                        20, 0, 20);

    SummaryStats st4;
    SummaryStats st8;
    fillHist(4, h4, st4);
    fillHist(8, h8, st8);

    if (st4.nEvents == 0 || st8.nEvents == 0)
    {
      std::cerr << "ERROR: no valid events found.\n";
      f->Close();
      return;
    }

    h4->SetLineColor(kBlue + 1);
    h8->SetLineColor(kRed + 1);
    h4->SetLineWidth(2);
    h8->SetLineWidth(2);
    h4->GetYaxis()->SetTitleOffset(2.0);
    h4->GetXaxis()->SetTitleOffset(1.1);

    double autoMax = std::max(h4->GetMaximum(), h8->GetMaximum());
    if (yMax > 0.0)
      h4->SetMaximum(yMax);
    else
      h4->SetMaximum(autoMax * 1.1);

    h4->Draw();
    h8->Draw("SAME");

    auto *leg = new TLegend(0.60, 0.78, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->AddEntry(h4, "4-neighbor", "l");
    leg->AddEntry(h8, "8-neighbor", "l");
    leg->Draw();

    auto *p = new TPaveText(0.55, 0.62, 0.88, 0.76, "NDC");
    p->SetFillColor(0);
    p->SetTextAlign(12);
    p->AddText(Form("thr: %.3f MeV", thrMeV));
    p->AddText(Form("cluster thr: %.3f MeV", thrClMeV));
    p->AddText(Form("Events: %lld", st4.nEvents));
    // p->Draw();
  }

  c->Update();

  // c->SaveAs(Form("nClusterSummary_%d.png", neighbor));
}

// ============================================================
// ① calarrのみ：nCluster==1 のとき
//    - total dE（225和）分布
//    - cluster size（その1クラスタのセグメント数）分布
//
// 引数:
//   fileName : ROOT file（TTree "calarr"）
//   nLayer,nSeg : 15,15（225）
//   thrMeV   : hit閾値（dE_MeV > thrMeV がクラスタ対象）
//   wrapSeg  : seg方向周期境界
//   neighbor : 4 or 8
//   totMaxMeV: total dE ヒストのx上限（適当に広めに）
// ============================================================
void drawOneClusterCal_TotAndSize(const char *fileName = "output/beta.root",
                                  int nLayer = 15,
                                  int nSeg = 15,
                                  double thrMeV = 1.0,
                                  bool wrapSeg = true,
                                  int neighbor = 8,
                                  double totMaxMeV = 200.0)
{
  if (neighbor != 4 && neighbor != 8)
  {
    std::cerr << "ERROR: neighbor must be 4 or 8 (given " << neighbor << ")\n";
    return;
  }

  const int nTot = nLayer * nSeg;

  TFile *f = TFile::Open(fileName, "READ");
  if (!f || f->IsZombie())
  {
    std::cerr << "ERROR: cannot open file: " << fileName << "\n";
    return;
  }

  TTree *t = (TTree *)f->Get("calarr");
  if (!t)
  {
    std::cerr << "ERROR: TTree 'calarr' not found.\n";
    f->Close();
    return;
  }

  int ev = -1;
  std::vector<double> *de = nullptr;
  t->SetBranchAddress("eventID", &ev);
  t->SetBranchAddress("dE_MeV", &de);

  const char *neighLabel = (neighbor == 4) ? "4-neigh" : "8-neigh";

  auto *hTot = new TH1D(Form("hCalTot_1cl_%d", neighbor),
                        Form("Cal total dE for nCluster=1 (%s, thr=%.3f MeV);Total dE (MeV);Events",
                             neighLabel, thrMeV),
                        200, 0.0, totMaxMeV);

  auto *hSize = new TH1I(Form("hClSize_1cl_%d", neighbor),
                         Form("Cluster size for nCluster=1 (%s, thr=%.3f MeV);Cluster size (#segments);Events",
                              neighLabel, thrMeV),
                         20, 0, 20);

  Long64_t nEnt = t->GetEntries();
  Long64_t nValid = 0;
  Long64_t nOne = 0;

  for (Long64_t i = 0; i < nEnt; ++i)
  {
    t->GetEntry(i);
    if (!de)
      continue;
    if ((int)de->size() < nTot)
      continue;
    ++nValid;

    // totalEdep (sum over 225, threshold関係なし)
    double totE = 0.0;
    for (int k = 0; k < nTot; ++k)
      totE += (*de)[k];

    // hit mask + visited
    std::vector<char> hit(nTot, 0), vis(nTot, 0);
    for (int k = 0; k < nTot; ++k)
    {
      if ((*de)[k] > thrMeV)
        hit[k] = 1;
    }

    int nCluster = 0;
    int firstClusterSize = 0;

    // BFS connected components
    for (int layer = 0; layer < nLayer; ++layer)
    {
      for (int seg = 0; seg < nSeg; ++seg)
      {
        const int start = idx2d(layer, seg, nSeg);
        if (!hit[start] || vis[start])
          continue;

        ++nCluster;

        int thisSize = 0;
        std::queue<std::pair<int, int>> q;
        q.push({layer, seg});
        vis[start] = 1;

        while (!q.empty())
        {
          auto [ly, sg] = q.front();
          q.pop();
          ++thisSize;

          for (int dly = -1; dly <= 1; ++dly)
          {
            const int nly = ly + dly;
            if (nly < 0 || nly >= nLayer)
              continue;

            for (int dsg = -1; dsg <= 1; ++dsg)
            {
              if (dly == 0 && dsg == 0)
                continue;
              if (neighbor == 4 && (dly != 0 && dsg != 0))
                continue;

              int nsg = sg + dsg;
              if (wrapSeg)
              {
                if (nsg >= nSeg)
                  nsg -= nSeg;
                if (nsg < 0)
                  nsg += nSeg;
              }
              else
              {
                if (nsg < 0 || nsg >= nSeg)
                  continue;
              }

              const int nb = idx2d(nly, nsg, nSeg);
              if (hit[nb] && !vis[nb])
              {
                vis[nb] = 1;
                q.push({nly, nsg});
              }
            }
          }
        }

        if (nCluster == 1)
          firstClusterSize = thisSize;
        if (nCluster > 1)
          break;
      }
      if (nCluster > 1)
        break;
    }

    if (nCluster == 1)
    {
      ++nOne;
      hTot->Fill(totE);
      hSize->Fill(firstClusterSize);
    }
  }

  if (nValid == 0)
  {
    std::cerr << "ERROR: no valid events found.\n";
    f->Close();
    return;
  }

  gStyle->SetOptStat(0);

  auto *c = new TCanvas(Form("cCalTotSize_1cl_%d", neighbor),
                        Form("calarr nCluster=1 (%s)", neighLabel),
                        900, 900);
  c->Divide(1, 2);

  c->cd(1);
  hTot->Draw();
  hTot->GetYaxis()->SetRangeUser(0,270);

  auto *p = new TPaveText(0.55, 0.62, 0.88, 0.88, "NDC");
  p->SetFillColor(0);
  p->SetTextAlign(12);
  p->AddText(Form("thr: %.3f MeV", thrMeV));
  p->AddText(Form("valid: %lld", nValid));
  p->AddText(Form("nCl=1: %lld (%.2f%%)", nOne, 100.0 * (double)nOne / (double)nValid));
  // p->Draw();

  c->cd(2);
  hSize->Draw();

  c->Update();
  // f は開いたまま（対話用）
}

// ============================================================
// ② targetのみ：nCluster==1 のときの Target dE 分布
//    - calarrで nCluster==1 を判定 → その eventID の target dE を詰める
//
// 前提（必要なら変更）:
//   calarr    : eventID(int), dE_MeV(vector<double>)
//   targethit : eventID(int), dE_MeV(double)
//
// 引数:
//   fileName  : ROOT file
//   nLayer,nSeg,thrMeV,wrapSeg,neighbor : クラスタ判定条件（calarrと同じ）
//   targetMaxMeV: target dE ヒストのx上限
// ============================================================
void drawOneClusterTargetE(const char *fileName = "output/beta.root",
                           int nLayer = 15,
                           int nSeg = 15,
                           double thrMeV = 1.0,
                           bool wrapSeg = true,
                           int neighbor = 4,
                           double targetMaxMeV = 100.0)
{
  if (neighbor != 4 && neighbor != 8)
  {
    std::cerr << "ERROR: neighbor must be 4 or 8 (given " << neighbor << ")\n";
    return;
  }

  const int nTot = nLayer * nSeg;

  TFile *f = TFile::Open(fileName, "READ");
  if (!f || f->IsZombie())
  {
    std::cerr << "ERROR: cannot open file: " << fileName << "\n";
    return;
  }

  TTree *tCal = (TTree *)f->Get("calarr");
  if (!tCal)
  {
    std::cerr << "ERROR: TTree 'calarr' not found.\n";
    f->Close();
    return;
  }

  // tree名が違うならここを変更
  TTree *tTgt = (TTree *)f->Get("target");
  if (!tTgt)
  {
    std::cerr << "ERROR: TTree 'target' not found.\n";
    std::cerr << "       (If your target tree name differs, change f->Get(\"target\")).\n";
    f->Close();
    return;
  }

  int evCal = -1;
  std::vector<double> *de = nullptr;
  tCal->SetBranchAddress("eventID", &evCal);
  tCal->SetBranchAddress("dE_MeV", &de);

  int evTgt = -1;
  double tgtE = 0.0;
  tTgt->SetBranchAddress("eventID", &evTgt);
  tTgt->SetBranchAddress("dE_MeV", &tgtE);

  // eventID -> target dE
  std::unordered_map<int, double> tgtByEv;
  tgtByEv.reserve((size_t)tTgt->GetEntries());

  const Long64_t nTgt = tTgt->GetEntries();
  for (Long64_t i = 0; i < nTgt; ++i)
  {
    tTgt->GetEntry(i);
    // もし同一eventが複数rowあり得るなら += に変える
    tgtByEv[evTgt] = tgtE;
  }

  const char *neighLabel = (neighbor == 4) ? "4-neigh" : "8-neigh";

  auto *hTgt = new TH1D(Form("hTargetE_1cl_%d", neighbor),
                        Form("Target dE for nCluster=1 (%s, thr=%.3f MeV);Target dE (MeV);Events",
                             neighLabel, thrMeV),
                        400, 0.0, targetMaxMeV);

  Long64_t nEnt = tCal->GetEntries();
  Long64_t nValid = 0;
  Long64_t nOne = 0;
  Long64_t nFill = 0;

  for (Long64_t i = 0; i < nEnt; ++i)
  {
    tCal->GetEntry(i);
    if (!de)
      continue;
    if ((int)de->size() < nTot)
      continue;
    ++nValid;

    // hit mask + visited
    std::vector<char> hit(nTot, 0), vis(nTot, 0);
    for (int k = 0; k < nTot; ++k)
    {
      if ((*de)[k] > thrMeV)
        hit[k] = 1;
    }

    int nCluster = 0;

    // BFS count only（sizeは不要）
    for (int layer = 0; layer < nLayer; ++layer)
    {
      for (int seg = 0; seg < nSeg; ++seg)
      {
        const int start = idx2d(layer, seg, nSeg);
        if (!hit[start] || vis[start])
          continue;

        ++nCluster;

        std::queue<std::pair<int, int>> q;
        q.push({layer, seg});
        vis[start] = 1;

        while (!q.empty())
        {
          auto [ly, sg] = q.front();
          q.pop();

          for (int dly = -1; dly <= 1; ++dly)
          {
            const int nly = ly + dly;
            if (nly < 0 || nly >= nLayer)
              continue;

            for (int dsg = -1; dsg <= 1; ++dsg)
            {
              if (dly == 0 && dsg == 0)
                continue;
              if (neighbor == 4 && (dly != 0 && dsg != 0))
                continue;

              int nsg = sg + dsg;
              if (wrapSeg)
              {
                if (nsg >= nSeg)
                  nsg -= nSeg;
                if (nsg < 0)
                  nsg += nSeg;
              }
              else
              {
                if (nsg < 0 || nsg >= nSeg)
                  continue;
              }

              const int nb = idx2d(nly, nsg, nSeg);
              if (hit[nb] && !vis[nb])
              {
                vis[nb] = 1;
                q.push({nly, nsg});
              }
            }
          }
        }

        if (nCluster > 1)
          break;
      }
      if (nCluster > 1)
        break;
    }

    if (nCluster == 1)
    {
      ++nOne;
      auto it = tgtByEv.find(evCal);
      if (it != tgtByEv.end())
      {
        hTgt->Fill(it->second);
        ++nFill;
      }
    }
  }

  if (nValid == 0)
  {
    std::cerr << "ERROR: no valid calarr events found.\n";
    f->Close();
    return;
  }

  gStyle->SetOptStat(0);

  auto *c = new TCanvas(Form("cTargetE_1cl_%d", neighbor),
                        Form("target dE (nCluster=1, %s)", neighLabel),
                        900, 650);

  hTgt->Draw();

  auto *p = new TPaveText(0.55, 0.62, 0.88, 0.88, "NDC");
  p->SetFillColor(0);
  p->SetTextAlign(12);
  p->AddText(Form("thr: %.3f MeV", thrMeV));
  p->AddText(Form("valid cal: %lld", nValid));
  p->AddText(Form("nCl=1: %lld (%.2f%%)", nOne, 100.0 * (double)nOne / (double)nValid));
  p->AddText(Form("target matched: %lld", nFill));
  // p->Draw();

  c->Update();
  // f は開いたまま（対話用）
}

// 1 event のクラスタを抽出（あなたの定義そのまま + クラスタ内idxリストも返す）
static void findClusters(const std::vector<double> &de,
                         int nLayer, int nSeg,
                         double thrHitMeV,
                         double thrClMeV,
                         bool wrapSeg,
                         int neighbor,
                         std::vector<std::vector<int>> &clusters,
                         std::vector<double> &clSumE)
{
  const int nTot = nLayer * nSeg;
  std::vector<char> hit(nTot, 0), vis(nTot, 0);
  for (int k = 0; k < nTot; ++k)
    if (de[k] > thrHitMeV)
      hit[k] = 1;

  clusters.clear();
  clSumE.clear();

  for (int layer = 0; layer < nLayer; ++layer)
  {
    for (int seg = 0; seg < nSeg; ++seg)
    {
      const int start = idx2d(layer, seg, nSeg);
      if (!hit[start] || vis[start])
        continue;

      std::vector<int> cl;
      double sumE = 0.0;

      std::queue<std::pair<int, int>> q;
      q.push({layer, seg});
      vis[start] = 1;

      while (!q.empty())
      {
        auto [ly, sg] = q.front();
        q.pop();
        const int idx = idx2d(ly, sg, nSeg);
        cl.push_back(idx);
        sumE += de[idx];

        for (int dly = -1; dly <= 1; ++dly)
        {
          const int nly = ly + dly;
          if (nly < 0 || nly >= nLayer)
            continue;
          for (int dsg = -1; dsg <= 1; ++dsg)
          {
            if (dly == 0 && dsg == 0)
              continue;
            if (neighbor == 4 && (dly != 0 && dsg != 0))
              continue;

            int nsg = sg + dsg;
            if (wrapSeg)
            {
              if (nsg >= nSeg)
                nsg -= nSeg;
              if (nsg < 0)
                nsg += nSeg;
            }
            else
            {
              if (nsg < 0 || nsg >= nSeg)
                continue;
            }

            const int nb = idx2d(nly, nsg, nSeg);
            if (hit[nb] && !vis[nb])
            {
              vis[nb] = 1;
              q.push({nly, nsg});
            }
          }
        }
      }

      if (thrClMeV >= 0.0 && sumE <= thrClMeV)
        continue;

      clusters.push_back(std::move(cl));
      clSumE.push_back(sumE);
    }
  }
}

void analyzeNonPiClusters(const char *fileName = "output/beta.root",
                          const char *hitTreeName = "calhit",
                          int nLayer = 15, int nSeg = 15,
                          double thrHitMeV = 1.0,
                          double thrClMeV = 1.0,
                          bool wrapSeg = true,
                          int neighbor = 8,
                          double sumEmax = 200.0,
                          int topN = 10)
{
  const int nTot = nLayer * nSeg;

  TFile *f = TFile::Open(fileName, "READ");
  if (!f || f->IsZombie())
  {
    std::cerr << "cannot open " << fileName << "\n";
    return;
  }

  TTree *tCal = (TTree *)f->Get("calarr");
  if (!tCal)
  {
    std::cerr << "no TTree calarr\n";
    return;
  }

  TTree *tHit = (TTree *)f->Get(hitTreeName);
  if (!tHit)
  {
    std::cerr << "no hit tree '" << hitTreeName << "'. まず TFile::ls() でツリー名確認して。\n";
    return;
  }

  // --- calarr
  int evCal = -1;
  std::vector<double> *de = nullptr;
  tCal->SetBranchAddress("eventID", &evCal);
  tCal->SetBranchAddress("dE_MeV", &de);

  // --- hit tree
  int evH = -1, copyNo = -1, pdg = 0;
  char creator[128] = {0};
  tHit->SetBranchAddress("eventID", &evH);
  tHit->SetBranchAddress("copyNo", &copyNo);
  tHit->SetBranchAddress("pdg", &pdg);
  tHit->SetBranchAddress("creator", creator); // char[128]

  // map (event,copy)->pdg/creator
  std::unordered_map<std::uint64_t, int> pdgByEvCopy;
  std::unordered_map<std::uint64_t, std::string> creatorByEvCopy;
  pdgByEvCopy.reserve((size_t)tHit->GetEntries() * 2);
  creatorByEvCopy.reserve((size_t)tHit->GetEntries() * 2);

  for (Long64_t i = 0; i < tHit->GetEntries(); ++i)
  {
    tHit->GetEntry(i);
    if (copyNo < 0 || copyNo >= nTot)
      continue;
    const auto key = keyEvCopy(evH, copyNo);
    pdgByEvCopy[key] = pdg;
    creatorByEvCopy[key] = std::string(creator);
  }

  // hist
  auto *hLeadCat = new TH1I("hLeadCat", "non-pi cluster lead category;category;clusters", 9, 0, 9);
  auto *hSumE_byCat = new TH2D("hSumE_byCat", "non-pi clusters;sumE (MeV);category",
                               300, 0, 100, 9, 0, 9);
  auto *hSumE_nonPi = new TH1D("hSumE_nonPi", "non-pi cluster sumE;sumE (MeV);clusters",
                               200, 0, sumEmax);
  auto *hSize_nonPi = new TH1I("hSize_nonPi", "non-pi cluster size;#segments;clusters",
                               50, 0, 50);

  // counters
  long long nCl_all = 0, nCl_pi = 0, nCl_nonPi = 0;
  std::unordered_map<int, long long> cntCat;
  std::unordered_map<int, double> sumEcat;

  // creator 集計：leadカテゴリごとに (creator -> count)
  std::unordered_map<int, std::unordered_map<std::string, long long>> cntCreatorByLead;
  struct ProcStat
  {
    long long nCl = 0;     // クラスタ数（クラスタ重み）
    double sumClE = 0.0;   // Σ cluster sumE
    double sumClE2 = 0.0;  // Σ cluster sumE^2（RMS用）
    long long sumSize = 0; // Σ cluster size

    // イベント重みで totE を平均したい場合用
    long long nEv = 0;       // そのcreatorが出たイベント数
    double sumTotE = 0.0;    // Σ event totE（同一eventで同一creatorは1回だけ）
    int lastEv = -999999999; // 直近に totE を足した eventID（tCal が event順ならこれで十分）
  };

  std::unordered_map<std::string, ProcStat> statByCreator;
  statByCreator.reserve(512);

  // loop
  std::vector<std::vector<int>> clusters;
  std::vector<double> clSumE;

  for (Long64_t iev = 0; iev < tCal->GetEntries(); ++iev)
  {
    tCal->GetEntry(iev);
    if (!de)
      continue;
    if ((int)de->size() < nTot)
      continue;

    double totE = 0.0;
    for (int k = 0; k < nTot; ++k)
      totE += (*de)[k];

    findClusters(*de, nLayer, nSeg, thrHitMeV, thrClMeV, wrapSeg, neighbor, clusters, clSumE);

    for (size_t ic = 0; ic < clusters.size(); ++ic)
    {
      ++nCl_all;

      const auto &cl = clusters[ic];
      const double sumE = clSumE[ic];

      bool hasPiMinus = false;
      double eByCat[9] = {0};

      // クラスタ代表creator = クラスタ内で dE 最大のセグメントの creator
      double bestSegE = -1.0;
      std::string repCreator = "Unknown";

      for (int idx : cl)
      {
        const int copy = idx; // 0..224 が copyNo と一致前提
        const auto key = keyEvCopy(evCal, copy);

        const int pdg_here = (pdgByEvCopy.count(key) ? pdgByEvCopy[key] : 0);
        if (pdg_here == -211)
          hasPiMinus = true;

        const int cat = pdgCategory(pdg_here);
        eByCat[cat] += (*de)[idx];

        const double segE = (*de)[idx];
        if (segE > bestSegE)
        {
          bestSegE = segE;
          auto itc = creatorByEvCopy.find(key);
          repCreator = (itc != creatorByEvCopy.end()) ? itc->second : "Unknown";
        }
      }

      if (hasPiMinus)
      {
        ++nCl_pi;
        continue;
      }

      ++nCl_nonPi;
      hSumE_nonPi->Fill(sumE);
      hSize_nonPi->Fill((int)cl.size());

      // --- per creator(process) stats ---
      auto &st = statByCreator[repCreator];
      st.nCl++;
      st.sumClE += sumE;
      st.sumClE2 += sumE * sumE;
      st.sumSize += (long long)cl.size();

      // イベント重みの平均 totE（同一イベント内で同一creatorは1回だけ totE を足す）
      if (st.lastEv != evCal)
      {
        st.lastEv = evCal;
        st.nEv++;
        st.sumTotE += totE;
      }

      // lead category = max energy in cluster
      int lead = kClCat_Other;
      double best = -1.0;
      for (int c = 1; c <= 8; ++c)
      {
        if (eByCat[c] > best)
        {
          best = eByCat[c];
          lead = c;
        }
      }

      // fill
      hLeadCat->Fill(lead);
      hSumE_byCat->Fill(sumE, lead);
      cntCat[lead] += 1;
      sumEcat[lead] += sumE;

      // creator count
      cntCreatorByLead[lead][repCreator] += 1;
    }
  }

  // summary
  std::cout << "clusters: all=" << nCl_all << "  pi-=" << nCl_pi << "  non-pi=" << nCl_nonPi << "\n";
  std::cout << "non-pi lead breakdown:\n";
  for (int c = 1; c <= 8; ++c)
  {
    const auto n = cntCat[c];
    const auto s = sumEcat[c];
    const double frac = (nCl_nonPi > 0) ? 100.0 * (double)n / (double)nCl_nonPi : 0.0;
    const double meanE = (n > 0) ? s / (double)n : 0.0;
    std::cout << "  " << catName(c)
              << " : N=" << n
              << " (" << frac << "%)"
              << "  meanSumE=" << meanE << " MeV\n";
  }

  auto printTop = [&](int lead)
  {
    auto it = cntCreatorByLead.find(lead);
    if (it == cntCreatorByLead.end() || it->second.empty())
    {
      std::cout << "\n[Lead " << catName(lead) << "] no creator entries\n";
      return;
    }
    auto &m = it->second;
    std::vector<std::pair<std::string, long long>> v(m.begin(), m.end());
    std::sort(v.begin(), v.end(),
              [](auto &a, auto &b)
              { return a.second > b.second; });

    std::cout << "\n[Lead " << catName(lead) << "] Top creators:\n";
    for (int i = 0; i < (int)v.size() && i < topN; ++i)
    {
      std::cout << "  " << v[i].first << " : " << v[i].second << "\n";
    }
  };

  // ---- print mean(cluster dE, size, event totE) per creator ----
  std::vector<std::pair<std::string, ProcStat>> vproc(statByCreator.begin(), statByCreator.end());
  std::sort(vproc.begin(), vproc.end(),
            [](auto &a, auto &b)
            { return a.second.nCl > b.second.nCl; });

  std::cout << "\n=== Per creator(process) mean dE ===\n";
  std::cout << "  (Cluster-weighted) meanClE = <sumE_cluster>\n";
  std::cout << "  (Event-weighted)   meanTotE = <sum over 225 in event> for events where creator appears\n\n";

  int nPrint = 0;
  for (auto &kv : vproc)
  {
    const auto &proc = kv.first;
    const auto &st = kv.second;
    if (st.nCl <= 0)
      continue;

    const double meanClE = st.sumClE / (double)st.nCl;
    const double varClE = st.sumClE2 / (double)st.nCl - meanClE * meanClE;
    const double rmsClE = (varClE > 0) ? std::sqrt(varClE) : 0.0;

    const double meanSize = (double)st.sumSize / (double)st.nCl;

    // event-weighted totE (if available)
    const double meanTotE = (st.nEv > 0) ? st.sumTotE / (double)st.nEv : 0.0;

    std::cout << "  " << proc
              << " : Ncl=" << st.nCl
              << "  meanClE=" << meanClE << " MeV"
              << "  rmsClE=" << rmsClE << " MeV"
              << "  meanSize=" << meanSize
              << "  Nev=" << st.nEv
              << "  meanTotE=" << meanTotE << " MeV\n";

    if (++nPrint >= topN)
      break;
  }

  for (int c = 1; c <= 8; ++c)
    printTop(c);

  // draw
  gStyle->SetOptStat(0);
  auto *c1 = new TCanvas("c_nonPi", "non-pi clusters", 1100, 800);
  c1->Divide(2, 2);
  c1->cd(1);
  hLeadCat->Draw();
  c1->cd(2);
  hSumE_nonPi->Draw();
  c1->cd(3);
  hSize_nonPi->Draw();
  c1->cd(4);
  hSumE_byCat->Draw("COLZ");
  gPad->SetLogz();

  for (int c = 1; c <= 8; ++c)
  {
    hLeadCat->GetXaxis()->SetBinLabel(c + 1, catName(c));
    hSumE_byCat->GetYaxis()->SetBinLabel(c + 1, catName(c));
  }
  c1->Update();
}
