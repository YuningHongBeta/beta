#include <TFile.h>
#include <TTree.h>

#include <cstdint>
#include <iostream>
#include <queue>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace {
int CountClusters(const std::vector<double> &de, int nLayer, int nSector,
                  double hitThreshold, double clusterThreshold)
{
  const int nCells = nLayer * nSector;
  if (static_cast<int>(de.size()) != nCells)
    return -1;
  std::vector<char> hit(nCells, 0), seen(nCells, 0);
  for (int i = 0; i < nCells; ++i)
    hit[i] = de[i] > hitThreshold;

  int nCluster = 0;
  for (int layer = 0; layer < nLayer; ++layer) {
    for (int sector = 0; sector < nSector; ++sector) {
      const int start = layer * nSector + sector;
      if (!hit[start] || seen[start])
        continue;
      double sum = 0.0;
      std::queue<std::pair<int, int>> pending;
      pending.push({layer, sector});
      seen[start] = 1;
      while (!pending.empty()) {
        const auto [currentLayer, currentSector] = pending.front();
        pending.pop();
        sum += de[currentLayer * nSector + currentSector];
        constexpr int dLayer[4] = {-1, 1, 0, 0};
        constexpr int dSector[4] = {0, 0, -1, 1};
        for (int direction = 0; direction < 4; ++direction) {
          const int nextLayer = currentLayer + dLayer[direction];
          if (nextLayer < 0 || nextLayer >= nLayer)
            continue;
          const int nextSector =
              (currentSector + dSector[direction] + nSector) % nSector;
          const int index = nextLayer * nSector + nextSector;
          if (hit[index] && !seen[index]) {
            seen[index] = 1;
            pending.push({nextLayer, nextSector});
          }
        }
      }
      if (sum > clusterThreshold)
        ++nCluster;
    }
  }
  return nCluster;
}
} // namespace

// Focused baseline summary. This deliberately retains the historical
// 4-neighbour rectangular clustering; it is not an optimized BGOegg classifier.
void summarize_bgo_pc(const char *path, double pcThresholdMeV = 0.5,
                      double hitThresholdMeV = 1.0,
                      double clusterThresholdMeV = 1.0)
{
  TFile file(path, "READ");
  auto *evt = dynamic_cast<TTree *>(file.Get("evt"));
  auto *calarr = dynamic_cast<TTree *>(file.Get("calarr"));
  auto *runmeta = dynamic_cast<TTree *>(file.Get("runmeta"));
  if (!evt || !calarr || !runmeta || !evt->GetBranch("EdepPC_MeV")) {
    std::cerr << "missing evt/calarr/runmeta or evt.EdepPC_MeV in " << path
              << '\n';
    return;
  }

  int nLayer = 0, nSector = 0;
  runmeta->SetBranchAddress("nLayer", &nLayer);
  runmeta->SetBranchAddress("nSector", &nSector);
  runmeta->GetEntry(0);

  int eventID = -1;
  double pcEdep = 0.0;
  evt->SetBranchAddress("eventID", &eventID);
  evt->SetBranchAddress("EdepPC_MeV", &pcEdep);
  std::unordered_map<int, double> pcByEvent;
  pcByEvent.reserve(evt->GetEntries());
  for (Long64_t entry = 0; entry < evt->GetEntries(); ++entry) {
    evt->GetEntry(entry);
    if (!pcByEvent.emplace(eventID, pcEdep).second)
      throw std::runtime_error("duplicate eventID in evt tree");
  }

  std::vector<double> *de = nullptr;
  calarr->SetBranchAddress("eventID", &eventID);
  calarr->SetBranchAddress("dE_MeV", &de);
  std::unordered_set<int> seenCalarr;
  seenCalarr.reserve(calarr->GetEntries());
  std::uint64_t nEvent = 0, nPCHit = 0, nOneCluster = 0, nOneClusterPCHit = 0;
  for (Long64_t entry = 0; entry < calarr->GetEntries(); ++entry) {
    calarr->GetEntry(entry);
    if (!seenCalarr.emplace(eventID).second)
      throw std::runtime_error("duplicate eventID in calarr tree");
    const auto found = pcByEvent.find(eventID);
    if (!de)
      throw std::runtime_error("null dE_MeV vector in calarr tree");
    if (found == pcByEvent.end())
      throw std::runtime_error("calarr eventID is missing from evt tree");
    const int nCluster = CountClusters(*de, nLayer, nSector,
                                       hitThresholdMeV,
                                       clusterThresholdMeV);
    if (nCluster < 0)
      throw std::runtime_error("calarr dE_MeV size disagrees with runmeta");
    ++nEvent;
    const bool pcHit = found->second > pcThresholdMeV;
    nPCHit += pcHit;
    if (nCluster == 1) {
      ++nOneCluster;
      nOneClusterPCHit += pcHit;
    }
  }
  if (seenCalarr.size() != pcByEvent.size())
    throw std::runtime_error("evt/calarr eventID sets differ");

  const double denominator = static_cast<double>(nEvent);
  const double oneClusterDenominator = static_cast<double>(nOneCluster);
  std::cout << "file=" << path
            << " cells=" << nLayer << 'x' << nSector
            << " events=" << nEvent
            << " pc_hit_pct=" << (denominator ? 100.0 * nPCHit / denominator : 0.0)
            << " ncluster1_pct=" << (denominator ? 100.0 * nOneCluster / denominator : 0.0)
            << " pc_veto_given_ncluster1_pct="
            << (oneClusterDenominator
                    ? 100.0 * nOneClusterPCHit / oneClusterDenominator
                    : 0.0)
            << " ncluster1_after_pc_pct="
            << (denominator
                    ? 100.0 * (nOneCluster - nOneClusterPCHit) / denominator
                    : 0.0)
            << '\n';
}
