#include <TBranch.h>
#include <TFile.h>
#include <TTree.h>

#include "../include/BGOeggGeometry.hh"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <queue>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr std::uint32_t kMagic = 0x42474f32;  // BGO2
constexpr std::uint32_t kVersion = 2;
constexpr int kNEfficiency = 4;
constexpr int kNPoissonFeaturePerEfficiency = 17;
constexpr double kEfficiencies[kNEfficiency] = {0.01, 0.02, 0.05, 0.10};
const char *kEfficiencyTags[kNEfficiency] = {"eff01", "eff02", "eff05", "eff10"};

enum Feature : int {
  kEventID = 0,
  kMaxE,
  kSumE,
  kNHit,
  kNCl4,
  kNCl8,
  kCl1Sum4,
  kCl2Sum4,
  kCl1Size4,
  kCl1MaxFrac4,
  kCl1RmsDeg4,
  kAllRmsDeg,
  kCl2OverCl1,
  kIsolatedEFrac4,
  kSeed15NCl4,
  kSeed15Cl1Sum4,
  kSeed15Cl1Size4,
  kE2OverE1,
  kNHit2T,
  kSum2T,
  kLocal3x3Frac,
  kMeanHitE,
  kThetaSpan,
  kPhiSpan,
  kTHPresent,
  kTHTotalE,
  kTHMaxE,
  kTHNHit,
  kTHLeadingSeg,
  kTHLeadingTimeNs,
  kTLCPresent,
  kTLCExpectedTotal,
  kTLCExpectedMax,
  kTLCNExpectedSeg,
  kTLCLeadingSeg,
  kTLCLeadingTimeNs,
  kMatchedExpectedPhotons,
  kMatchedExpectedFraction,
  kMatchedTLCMaxExpected,
  kMatchedTLCLeadingSeg,
  kCherenkovDtNs,
  kCherenkovDtValid,
  kBGOLeadingCell,
  kBGOLeadingThetaDeg,
  kBGOLeadingPhiDeg,
  kFinalReadoutPresent,
  kTHBgoMatchedSeg,
  kTHBgoMatchedE,
  kTHGeomPathMm,
  kTHMatchedDEdx,
  kTHMatchedZRecoMm,
  kTHMatchedZRecoValid,
  kZPredMm,
  kAbsDeltaZMm,
  kAbsDeltaZValid,
  kAbsDeltaZLt90,
  kAbsDeltaZSigT02,
  kAbsDeltaZLt90SigT02,
  kAbsDeltaZSigT05,
  kAbsDeltaZLt90SigT05,
  kAbsDeltaZSigT10,
  kAbsDeltaZLt90SigT10,
  kPCSumE,
  kPCDownE,
  kPCUpE,
  kPoissonFeatureStart,
  kNFeature =
      kPoissonFeatureStart + kNEfficiency * kNPoissonFeaturePerEfficiency
};

const std::vector<std::string> kFeatureNames = [] {
  std::vector<std::string> names = {
      "eventID", "maxE", "sumE", "nHit", "nCl4", "nCl8", "cl1Sum4",
      "cl2Sum4", "cl1Size4", "cl1MaxFrac4", "cl1RmsDeg4", "allRmsDeg",
      "cl2OverCl1", "isolatedEFrac4", "seed15NCl4", "seed15Cl1Sum4",
      "seed15Cl1Size4", "e2OverE1", "nHit2T", "sum2T", "local3x3Frac",
      "meanHitE", "thetaSpan", "phiSpan", "thPresent", "thTotalE_MeV",
      "thMaxE_MeV", "thNHit", "thLeadingSeg", "thLeadingTime_ns",
      "tlcPresent", "tlcExpectedTotal", "tlcExpectedMax", "tlcNExpectedSeg",
      "tlcLeadingSeg", "tlcLeadingTime_ns", "matchedExpectedPhotons",
      "matchedExpectedFraction", "matchedTlcMaxExpected", "matchedTlcLeadingSeg",
      "cherenkovDt_ns", "cherenkovDtValid", "bgoLeadingCell",
      "bgoLeadingTheta_deg", "bgoLeadingPhi_deg", "finalReadoutPresent",
      "thBgoMatchedSeg", "thBgoMatchedE_MeV", "thGeomPath_mm",
      "thMatchedDEdx_MeV_per_mm", "thMatchedZReco_mm",
      "thMatchedZRecoValid", "zPred_mm", "absDeltaZ_mm",
      "absDeltaZValid", "absDeltaZLt90", "absDeltaZ_sigT0p2ns_mm",
      "absDeltaZLt90_sigT0p2ns", "absDeltaZ_sigT0p5ns_mm",
      "absDeltaZLt90_sigT0p5ns", "absDeltaZ_sigT1p0ns_mm",
      "absDeltaZLt90_sigT1p0ns", "pcSumE_MeV", "pcDownE_MeV",
      "pcUpE_MeV"};
  for (int ieff = 0; ieff < kNEfficiency; ++ieff) {
    const std::string tag(kEfficiencyTags[ieff]);
    names.push_back("tlcNpeTotal_" + tag);
    names.push_back("tlcNpeMax_" + tag);
    names.push_back("tlcNpeNSegGe1_" + tag);
    names.push_back("tlcNpeNSegGe2_" + tag);
    names.push_back("tlcNpeNSegGe3_" + tag);
    names.push_back("matchedNpeTotal_" + tag);
    names.push_back("matchedNpeFraction_" + tag);
    names.push_back("matchedNpeMax_" + tag);
    names.push_back("matchedNpeNSegGe1_" + tag);
    names.push_back("matchedNpeNSegGe2_" + tag);
    names.push_back("matchedNpeNSegGe3_" + tag);
    names.push_back("tlcAnySegHitGe1_" + tag);
    names.push_back("tlcAnySegHitGe2_" + tag);
    names.push_back("tlcAnySegHitGe3_" + tag);
    names.push_back("matchedAnySegHitGe1_" + tag);
    names.push_back("matchedAnySegHitGe2_" + tag);
    names.push_back("matchedAnySegHitGe3_" + tag);
  }
  return names;
}();

int poissonFeature(int efficiencyIndex, int offset) {
  return kPoissonFeatureStart +
         efficiencyIndex * kNPoissonFeaturePerEfficiency + offset;
}

std::uint64_t splitmix64(std::uint64_t value) {
  value += 0x9e3779b97f4a7c15ULL;
  value = (value ^ (value >> 30U)) * 0xbf58476d1ce4e5b9ULL;
  value = (value ^ (value >> 27U)) * 0x94d049bb133111ebULL;
  return value ^ (value >> 31U);
}

struct Meta {
  int nLayer = 0;
  int nSector = 0;
  int segmentationMode = -1;
  int geometryMode = -1;
  int physicsFlag = -1;
  int nSegTH = 0;
  int nSegTLC = 0;
  double thetaMinDeg = 0.0;
  double thetaMaxDeg = 0.0;
  double rMinCm = 0.0;
  double thicknessCm = 0.0;
  double bgoZOffsetCm = 0.0;
  double thRMinMm = 15.0;
  double thRMaxMm = 21.0;
  double thBarZMinMm = -300.0;
  double thBarZMaxMm = 300.0;
  double thEffectiveLightSpeedMmPerNs = 150.0;
  int thTimingSmearingApplied = 0;
  bool finalReadoutMetadata = false;
};

struct Geometry {
  Meta meta;
  std::vector<std::array<double, 3>> directions;

  explicit Geometry(const Meta &m) : meta(m) {
    if (meta.nLayer <= 0 || meta.nSector <= 0)
      throw std::runtime_error("invalid nLayer/nSector");
    if (meta.segmentationMode < 0 || meta.segmentationMode > 2)
      throw std::runtime_error("unknown segmentationMode");
    directions.reserve(NCell());

    if (meta.segmentationMode == 2) {
      if (meta.geometryMode != 2)
        throw std::runtime_error(
            "bgoegg_published segmentation requires geometryMode=2");
      const auto rings = BGOeggGeometry::BuildRings(meta.nLayer);
      if (static_cast<int>(rings.size()) != meta.nLayer || meta.nSector != 60)
        throw std::runtime_error("invalid published BGOegg ring/sector count");
      const double phiHalf =
          0.5 * BGOeggGeometry::kPhiSpanDeg * kPi / 180.0;
      for (const auto &ring : rings) {
        const double thetaLow = ring.thetaLowDeg * kPi / 180.0;
        const double thetaHigh = ring.thetaHighDeg * kPi / 180.0;
        const double front = ring.frontRadiusMm;
        const double rear = front + BGOeggGeometry::kCrystalLengthMm;
        const double radialSum = front + rear;
        const double rhoCenter = 0.25 * radialSum *
            (std::sin(thetaLow) + std::sin(thetaHigh)) * std::cos(phiHalf);
        const double zCenterLocal = 0.25 * radialSum *
            (std::cos(thetaLow) + std::cos(thetaHigh));
        const double zCenter = zCenterLocal + 10.0 * meta.bgoZOffsetCm;
        for (int ip = 0; ip < meta.nSector; ++ip) {
          const double phi = 2.0 * kPi * (ip + 0.5) / meta.nSector;
          const double x = rhoCenter * std::cos(phi);
          const double y = rhoCenter * std::sin(phi);
          const double norm = std::sqrt(x * x + y * y + zCenter * zCenter);
          if (!(norm > 0.0))
            throw std::runtime_error("invalid BGOegg cell centroid direction");
          directions.push_back({x / norm, y / norm, zCenter / norm});
        }
      }
      return;
    }

    const double thetaMin = meta.thetaMinDeg * kPi / 180.0;
    const double thetaMax = meta.thetaMaxDeg * kPi / 180.0;
    const double cosMin = std::cos(thetaMin);
    const double cosMax = std::cos(thetaMax);
    for (int it = 0; it < meta.nLayer; ++it) {
      double thetaCenter = 0.0;
      if (meta.segmentationMode == 1) {
        const double c0 = cosMin + (cosMax - cosMin) *
            (static_cast<double>(it) / meta.nLayer);
        const double c1 = cosMin + (cosMax - cosMin) *
            (static_cast<double>(it + 1) / meta.nLayer);
        thetaCenter = std::acos(std::clamp(0.5 * (c0 + c1), -1.0, 1.0));
      } else {
        const double t0 = thetaMin + (thetaMax - thetaMin) *
            (static_cast<double>(it) / meta.nLayer);
        const double t1 = thetaMin + (thetaMax - thetaMin) *
            (static_cast<double>(it + 1) / meta.nLayer);
        thetaCenter = 0.5 * (t0 + t1);
      }
      for (int ip = 0; ip < meta.nSector; ++ip) {
        const double phi = 2.0 * kPi * (ip + 0.5) / meta.nSector;
        const double radiusMm =
            10.0 * (meta.rMinCm + 0.5 * meta.thicknessCm);
        const double x = radiusMm * std::sin(thetaCenter) * std::cos(phi);
        const double y = radiusMm * std::sin(thetaCenter) * std::sin(phi);
        const double z = radiusMm * std::cos(thetaCenter) +
                         10.0 * meta.bgoZOffsetCm;
        const double norm = std::sqrt(x * x + y * y + z * z);
        directions.push_back({x / norm, y / norm, z / norm});
      }
    }
  }

  int NCell() const { return meta.nLayer * meta.nSector; }
  int Index(int theta, int phi) const { return theta * meta.nSector + phi; }
};

struct Cluster {
  std::vector<int> cells;
  double sum = 0.0;
  double max = 0.0;
  bool hasSeed15 = false;
};

void requireBranch(TTree *tree, const char *name) {
  if (!tree || !tree->GetBranch(name))
    throw std::runtime_error(std::string("missing branch ") +
                             (tree ? tree->GetName() : "<null>") + "." + name);
}

Meta readMeta(TTree *tree, bool allowLegacyHodo) {
  if (!tree || tree->GetEntries() != 1)
    throw std::runtime_error("runmeta must contain exactly one entry");
  const char *required[] = {"nLayer", "nSector", "segmentationMode",
                            "geometryMode", "physicsFlag", "thetaMin_deg",
                            "thetaMax_deg", "rMin_cm", "thickness_cm",
                            "nSegTH", "nSegTLC"};
  for (const char *name : required) requireBranch(tree, name);
  const char *finalRequired[] = {"thRMin_mm", "thRMax_mm", "thBarZMin_mm",
                                 "thBarZMax_mm",
                                 "thEffectiveLightSpeed_mm_per_ns"};
  bool finalMetadata = true;
  for (const char *name : finalRequired)
    finalMetadata = finalMetadata && tree->GetBranch(name);
  const char *finalMetadataPresence[] = {"thTimingSmearingApplied",
                                         "thTimingModel"};
  for (const char *name : finalMetadataPresence)
    finalMetadata = finalMetadata && tree->GetBranch(name);
  if (!finalMetadata && !allowLegacyHodo)
    throw std::runtime_error(
        "final TH runmeta missing; use --allow-legacy-hodo only for old exploratory files");
  Meta m;
  m.finalReadoutMetadata = finalMetadata;
  tree->SetBranchStatus("*", 0);
  for (const char *name : required) tree->SetBranchStatus(name, 1);
  if (finalMetadata)
    for (const char *name : finalRequired) tree->SetBranchStatus(name, 1);
  if (finalMetadata) tree->SetBranchStatus("thTimingSmearingApplied", 1);
  tree->SetBranchAddress("nLayer", &m.nLayer);
  tree->SetBranchAddress("nSector", &m.nSector);
  tree->SetBranchAddress("segmentationMode", &m.segmentationMode);
  tree->SetBranchAddress("geometryMode", &m.geometryMode);
  tree->SetBranchAddress("physicsFlag", &m.physicsFlag);
  tree->SetBranchAddress("thetaMin_deg", &m.thetaMinDeg);
  tree->SetBranchAddress("thetaMax_deg", &m.thetaMaxDeg);
  tree->SetBranchAddress("rMin_cm", &m.rMinCm);
  tree->SetBranchAddress("thickness_cm", &m.thicknessCm);
  tree->SetBranchAddress("nSegTH", &m.nSegTH);
  tree->SetBranchAddress("nSegTLC", &m.nSegTLC);
  if (tree->GetBranch("bgoZOffset_cm")) {
    tree->SetBranchStatus("bgoZOffset_cm", 1);
    tree->SetBranchAddress("bgoZOffset_cm", &m.bgoZOffsetCm);
  }
  if (finalMetadata) {
    tree->SetBranchAddress("thRMin_mm", &m.thRMinMm);
    tree->SetBranchAddress("thRMax_mm", &m.thRMaxMm);
    tree->SetBranchAddress("thBarZMin_mm", &m.thBarZMinMm);
    tree->SetBranchAddress("thBarZMax_mm", &m.thBarZMaxMm);
    tree->SetBranchAddress("thEffectiveLightSpeed_mm_per_ns",
                           &m.thEffectiveLightSpeedMmPerNs);
    tree->SetBranchAddress("thTimingSmearingApplied",
                           &m.thTimingSmearingApplied);
  }
  tree->GetEntry(0);
  if (!(m.thetaMinDeg < m.thetaMaxDeg) || !(m.rMinCm > 0.0) ||
      !(m.thicknessCm > 0.0) || m.nSegTH <= 0 || m.nSegTLC <= 0)
    throw std::runtime_error("invalid runmeta geometry/hodoscope values");
  if (!(m.thRMinMm > 0.0 && m.thRMinMm < m.thRMaxMm &&
        m.thBarZMinMm < m.thBarZMaxMm &&
        m.thEffectiveLightSpeedMmPerNs > 0.0))
    throw std::runtime_error("invalid final TH geometry/timing metadata");
  if (m.finalReadoutMetadata && m.thTimingSmearingApplied != 0)
    throw std::runtime_error(
        "strict timing response scan requires unsmeared source TH times");
  return m;
}

std::unordered_map<int, Long64_t> buildEventIndex(TTree *tree) {
  requireBranch(tree, "eventID");
  tree->SetBranchStatus("*", 0);
  tree->SetBranchStatus("eventID", 1);
  int eventID = -1;
  tree->SetBranchAddress("eventID", &eventID);
  std::unordered_map<int, Long64_t> index;
  index.reserve(static_cast<size_t>(tree->GetEntries() * 1.3 + 1));
  for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
    tree->GetEntry(i);
    if (!index.emplace(eventID, i).second)
      throw std::runtime_error(std::string("duplicate eventID in ") +
                               tree->GetName() + ": " + std::to_string(eventID));
  }
  return index;
}

std::vector<Cluster> findClusters(const Geometry &g,
                                  const std::vector<double> &e,
                                  double threshold, bool diagonal) {
  std::vector<unsigned char> visited(g.NCell(), 0);
  std::vector<Cluster> out;
  for (int start = 0; start < g.NCell(); ++start) {
    if (visited[start] || e[start] < threshold) continue;
    Cluster cl;
    std::queue<int> q;
    q.push(start);
    visited[start] = 1;
    while (!q.empty()) {
      const int idx = q.front();
      q.pop();
      cl.cells.push_back(idx);
      cl.sum += e[idx];
      cl.max = std::max(cl.max, e[idx]);
      if (e[idx] >= 1.5 * threshold) cl.hasSeed15 = true;
      const int it = idx / g.meta.nSector;
      const int ip = idx % g.meta.nSector;
      for (int dt = -1; dt <= 1; ++dt) {
        const int nt = it + dt;
        if (nt < 0 || nt >= g.meta.nLayer) continue;
        for (int dp = -1; dp <= 1; ++dp) {
          if (dt == 0 && dp == 0) continue;
          if (!diagonal && dt != 0 && dp != 0) continue;
          const int np = (ip + dp + g.meta.nSector) % g.meta.nSector;
          const int ni = g.Index(nt, np);
          if (!visited[ni] && e[ni] >= threshold) {
            visited[ni] = 1;
            q.push(ni);
          }
        }
      }
    }
    out.push_back(std::move(cl));
  }
  std::sort(out.begin(), out.end(), [](const Cluster &a, const Cluster &b) {
    return a.sum > b.sum;
  });
  return out;
}

double angularRmsDeg(const Geometry &g, const std::vector<double> &e,
                     const std::vector<int> &cells) {
  double sx = 0.0, sy = 0.0, sz = 0.0, sum = 0.0;
  for (int idx : cells) {
    const auto &u = g.directions[idx];
    sx += e[idx] * u[0];
    sy += e[idx] * u[1];
    sz += e[idx] * u[2];
    sum += e[idx];
  }
  const double norm = std::sqrt(sx * sx + sy * sy + sz * sz);
  if (sum <= 0.0 || norm <= 0.0) return 0.0;
  sx /= norm;
  sy /= norm;
  sz /= norm;
  double var = 0.0;
  for (int idx : cells) {
    const auto &u = g.directions[idx];
    const double dot =
        std::clamp(sx * u[0] + sy * u[1] + sz * u[2], -1.0, 1.0);
    const double angle = std::acos(dot) * 180.0 / kPi;
    var += e[idx] * angle * angle;
  }
  return std::sqrt(var / sum);
}

std::vector<float> makeBgoFeatures(const Geometry &g,
                                   const std::vector<double> &raw,
                                   int eventID, double threshold) {
  std::vector<double> e(g.NCell(), 0.0);
  std::vector<int> hitCells;
  std::vector<double> sortedE;
  double sum = 0.0, maxE = 0.0, sum2t = 0.0;
  int nHit2t = 0, thetaMin = g.meta.nLayer, thetaMax = -1;
  int leadingCell = -1;
  for (int i = 0; i < g.NCell(); ++i) {
    if (raw[i] < threshold) continue;
    e[i] = raw[i];
    hitCells.push_back(i);
    sortedE.push_back(e[i]);
    sum += e[i];
    if (e[i] > maxE) {
      maxE = e[i];
      leadingCell = i;
    }
    if (e[i] >= 2.0 * threshold) {
      ++nHit2t;
      sum2t += e[i];
    }
    const int it = i / g.meta.nSector;
    thetaMin = std::min(thetaMin, it);
    thetaMax = std::max(thetaMax, it);
  }
  std::sort(sortedE.begin(), sortedE.end(), std::greater<double>());
  auto cl4 = findClusters(g, e, threshold, false);
  auto cl8 = findClusters(g, e, threshold, true);
  std::vector<Cluster> seeded;
  double isolated = 0.0;
  for (const auto &cl : cl4) {
    if (cl.cells.size() == 1) isolated += cl.sum;
    if (cl.hasSeed15) seeded.push_back(cl);
  }
  const double cl1 = cl4.empty() ? 0.0 : cl4[0].sum;
  const double cl2 = cl4.size() < 2 ? 0.0 : cl4[1].sum;
  double local = 0.0;
  if (!hitCells.empty()) {
    const int imax =
        static_cast<int>(std::max_element(e.begin(), e.end()) - e.begin());
    const int mt = imax / g.meta.nSector;
    const int mp = imax % g.meta.nSector;
    for (int dt = -1; dt <= 1; ++dt) {
      if (mt + dt < 0 || mt + dt >= g.meta.nLayer) continue;
      for (int dp = -1; dp <= 1; ++dp) {
        const int pp = (mp + dp + g.meta.nSector) % g.meta.nSector;
        local += e[g.Index(mt + dt, pp)];
      }
    }
  }
  int phiSpan = 0;
  if (!hitCells.empty()) {
    std::vector<unsigned char> phiHit(g.meta.nSector, 0);
    for (int idx : hitCells) phiHit[idx % g.meta.nSector] = 1;
    int longestEmpty = 0, currentEmpty = 0;
    for (int k = 0; k < 2 * g.meta.nSector; ++k) {
      if (!phiHit[k % g.meta.nSector])
        currentEmpty = std::min(currentEmpty + 1, g.meta.nSector);
      else
        currentEmpty = 0;
      longestEmpty = std::max(longestEmpty, currentEmpty);
    }
    phiSpan = g.meta.nSector - longestEmpty;
  }
  std::vector<float> f(kNFeature, 0.0f);
  f[kEventID] = eventID;
  f[kMaxE] = maxE;
  f[kSumE] = sum;
  f[kNHit] = hitCells.size();
  f[kNCl4] = cl4.size();
  f[kNCl8] = cl8.size();
  f[kCl1Sum4] = cl1;
  f[kCl2Sum4] = cl2;
  f[kCl1Size4] = cl4.empty() ? 0.0 : cl4[0].cells.size();
  f[kCl1MaxFrac4] = cl1 > 0.0 ? cl4[0].max / cl1 : 0.0;
  f[kCl1RmsDeg4] =
      cl4.empty() ? 0.0 : angularRmsDeg(g, e, cl4[0].cells);
  f[kAllRmsDeg] = angularRmsDeg(g, e, hitCells);
  f[kCl2OverCl1] = cl1 > 0.0 ? cl2 / cl1 : 0.0;
  f[kIsolatedEFrac4] = sum > 0.0 ? isolated / sum : 0.0;
  f[kSeed15NCl4] = seeded.size();
  f[kSeed15Cl1Sum4] = seeded.empty() ? 0.0 : seeded[0].sum;
  f[kSeed15Cl1Size4] = seeded.empty() ? 0.0 : seeded[0].cells.size();
  f[kE2OverE1] =
      sortedE.size() >= 2 && sortedE[0] > 0.0 ? sortedE[1] / sortedE[0] : 0.0;
  f[kNHit2T] = nHit2t;
  f[kSum2T] = sum2t;
  f[kLocal3x3Frac] = sum > 0.0 ? local / sum : 0.0;
  f[kMeanHitE] = hitCells.empty() ? 0.0 : sum / hitCells.size();
  f[kThetaSpan] = hitCells.empty() ? 0.0 : thetaMax - thetaMin + 1;
  f[kPhiSpan] = phiSpan;
  f[kTHLeadingSeg] = -1.0f;
  f[kTHLeadingTimeNs] = -1.0f;
  f[kTLCLeadingSeg] = -1.0f;
  f[kTLCLeadingTimeNs] = -1.0f;
  f[kMatchedTLCLeadingSeg] = -1.0f;
  f[kBGOLeadingCell] = leadingCell;
  f[kTHBgoMatchedSeg] = -1.0f;
  f[kTHMatchedZRecoMm] = 0.0f;
  if (leadingCell >= 0) {
    const auto &u = g.directions[leadingCell];
    f[kBGOLeadingThetaDeg] =
        std::acos(std::clamp(u[2], -1.0, 1.0)) * 180.0 / kPi;
    double phi = std::atan2(u[1], u[0]) * 180.0 / kPi;
    if (phi < 0.0) phi += 360.0;
    f[kBGOLeadingPhiDeg] = phi;
  }
  return f;
}

void addHodoscopeFeatures(std::vector<float> &f, const Meta &m,
                          int eventID,
                          const std::vector<double> *thE,
                          const std::vector<double> *thTime,
                          const std::vector<double> *thZReco,
                          const std::vector<double> *tlcExpected,
                          const std::vector<double> *tlcTime) {
  if (thE && thTime) {
    f[kTHPresent] = 1.0f;
    double total = 0.0, maxE = 0.0;
    int nHit = 0, leading = -1;
    for (int i = 0; i < m.nSegTH; ++i) {
      total += (*thE)[i];
      if ((*thE)[i] > 0.0) ++nHit;
      if ((*thE)[i] > maxE) {
        maxE = (*thE)[i];
        leading = i;
      }
    }
    f[kTHTotalE] = total;
    f[kTHMaxE] = maxE;
    f[kTHNHit] = nHit;
    f[kTHLeadingSeg] = leading;
    if (leading >= 0) f[kTHLeadingTimeNs] = (*thTime)[leading];
  }
  if (tlcExpected && tlcTime) {
    f[kTLCPresent] = 1.0f;
    double total = 0.0, maxExpected = 0.0;
    int nSeg = 0, leading = -1;
    for (int i = 0; i < m.nSegTLC; ++i) {
      total += (*tlcExpected)[i];
      if ((*tlcExpected)[i] > 0.0) ++nSeg;
      if ((*tlcExpected)[i] > maxExpected) {
        maxExpected = (*tlcExpected)[i];
        leading = i;
      }
    }
    f[kTLCExpectedTotal] = total;
    f[kTLCExpectedMax] = maxExpected;
    f[kTLCNExpectedSeg] = nSeg;
    f[kTLCLeadingSeg] = leading;
    if (leading >= 0) f[kTLCLeadingTimeNs] = (*tlcTime)[leading];
  }

  const bool finalReadout =
      m.finalReadoutMetadata && thE && thTime && thZReco;
  f[kFinalReadoutPresent] = finalReadout ? 1.0f : 0.0f;

  // Extrapolate the target-center -> leading-BGO-cell line to the TH shell.
  // The TH segment used for dE/dx and z matching is the largest deposit among
  // the predicted azimuthal segment and its immediate neighbours.
  int thBgoMatched = -1;
  const int bgoLeading = static_cast<int>(f[kBGOLeadingCell]);
  if (finalReadout && bgoLeading >= 0) {
    double phi = f[kBGOLeadingPhiDeg] * kPi / 180.0;
    int predicted = static_cast<int>(std::floor(
        phi / (2.0 * kPi / static_cast<double>(m.nSegTH))));
    predicted = std::clamp(predicted, 0, m.nSegTH - 1);
    double matchedE = -1.0;
    for (int delta = -1; delta <= 1; ++delta) {
      const int seg = (predicted + delta + m.nSegTH) % m.nSegTH;
      if ((*thE)[seg] > matchedE) {
        matchedE = (*thE)[seg];
        thBgoMatched = seg;
      }
    }
    f[kTHBgoMatchedSeg] = thBgoMatched;
    f[kTHBgoMatchedE] = std::max(0.0, matchedE);

    const double theta = f[kBGOLeadingThetaDeg] * kPi / 180.0;
    const double sinTheta = std::sin(theta);
    const double cosTheta = std::cos(theta);
    if (sinTheta > 1.0e-12) {
      const double cotTheta = cosTheta / sinTheta;
      const double zInner = m.thRMinMm * cotTheta;
      const double zOuter = m.thRMaxMm * cotTheta;
      const bool crossesBarrel =
          zInner >= m.thBarZMinMm && zInner <= m.thBarZMaxMm &&
          zOuter >= m.thBarZMinMm && zOuter <= m.thBarZMaxMm;
      if (crossesBarrel) {
        const double path = (m.thRMaxMm - m.thRMinMm) / sinTheta;
        const double zPred = 0.5 * (m.thRMinMm + m.thRMaxMm) * cotTheta;
        f[kTHGeomPathMm] = path;
        f[kTHMatchedDEdx] = path > 0.0 ? f[kTHBgoMatchedE] / path : 0.0;
        f[kZPredMm] = zPred;
        if (thBgoMatched >= 0 && std::isfinite((*thZReco)[thBgoMatched]) &&
            (*thZReco)[thBgoMatched] >= m.thBarZMinMm &&
            (*thZReco)[thBgoMatched] <= m.thBarZMaxMm &&
            (*thE)[thBgoMatched] > 0.0) {
          const double zReco = (*thZReco)[thBgoMatched];
          const double absDeltaZ = std::abs(zReco - zPred);
          f[kTHMatchedZRecoMm] = zReco;
          f[kTHMatchedZRecoValid] = 1.0f;
          f[kAbsDeltaZMm] = absDeltaZ;
          f[kAbsDeltaZValid] = 1.0f;
          f[kAbsDeltaZLt90] = absDeltaZ < 90.0 ? 1.0f : 0.0f;

          // Provisional end-time resolution scan. Each end gets an independent
          // Gaussian fluctuation; the event/segment/sigma seed has no sample
          // label. Source timing must be unsmeared (enforced in readMeta).
          constexpr double sigmasNs[] = {0.2, 0.5, 1.0};
          constexpr int absFeatures[] = {
              kAbsDeltaZSigT02, kAbsDeltaZSigT05, kAbsDeltaZSigT10};
          constexpr int cutFeatures[] = {
              kAbsDeltaZLt90SigT02, kAbsDeltaZLt90SigT05,
              kAbsDeltaZLt90SigT10};
          for (int isigma = 0; isigma < 3; ++isigma) {
            const std::uint64_t seed = splitmix64(
                0x54485a534d454152ULL ^
                (static_cast<std::uint64_t>(
                     static_cast<std::uint32_t>(eventID)) << 16U) ^
                (static_cast<std::uint64_t>(thBgoMatched) << 4U) ^
                static_cast<std::uint64_t>(isigma));
            std::mt19937_64 random(seed);
            std::normal_distribution<double> gaussian(0.0, sigmasNs[isigma]);
            const double deltaLeft = gaussian(random);
            const double deltaRight = gaussian(random);
            const double smearedZ = zReco +
                0.5 * m.thEffectiveLightSpeedMmPerNs *
                (deltaLeft - deltaRight);
            const double smearedAbsDeltaZ = std::abs(smearedZ - zPred);
            f[absFeatures[isigma]] = smearedAbsDeltaZ;
            f[cutFeatures[isigma]] =
                smearedAbsDeltaZ < 90.0 ? 1.0f : 0.0f;
          }
        }
      }
    }
  }

  const int thLeading = thBgoMatched >= 0
                            ? thBgoMatched
                            : static_cast<int>(f[kTHLeadingSeg]);
  if (!tlcExpected || !tlcTime) return;

  // Match TLC segment centers to the azimuthal window covered by the leading
  // BGO-matched TH segment and its two neighbors. For equal segmentation this is exactly
  // {leading-1, leading, leading+1}; it also remains defined if counts differ.
  std::vector<unsigned char> matchedMask(m.nSegTLC, 0);
  double matched = 0.0, matchedMax = 0.0;
  int matchedLeading = -1;
  if (thLeading >= 0) {
    const double thCenter = 2.0 * kPi * (thLeading + 0.5) / m.nSegTH;
    const double halfWindow = 1.5 * (2.0 * kPi / m.nSegTH);
    for (int i = 0; i < m.nSegTLC; ++i) {
      const double center = 2.0 * kPi * (i + 0.5) / m.nSegTLC;
      const double delta = std::remainder(center - thCenter, 2.0 * kPi);
      if (std::abs(delta) > halfWindow + 1e-12) continue;
      matchedMask[i] = 1;
      matched += (*tlcExpected)[i];
      if ((*tlcExpected)[i] > matchedMax) {
        matchedMax = (*tlcExpected)[i];
        matchedLeading = i;
      }
    }
  }
  f[kMatchedExpectedPhotons] = matched;
  f[kMatchedExpectedFraction] =
      f[kTLCExpectedTotal] > 0.0f ? matched / f[kTLCExpectedTotal] : 0.0;
  f[kMatchedTLCMaxExpected] = matchedMax;
  f[kMatchedTLCLeadingSeg] = matchedLeading;
  if (matchedLeading >= 0 && f[kTHLeadingTimeNs] >= 0.0f &&
      (*tlcTime)[matchedLeading] >= 0.0) {
    f[kCherenkovDtNs] = (*tlcTime)[matchedLeading] - f[kTHLeadingTimeNs];
    f[kCherenkovDtValid] = 1.0f;
  }

  // Detector-response scan: turn analytic Frank-Tamm expected photons into
  // observed photoelectrons with a reproducible event/efficiency-specific
  // Poisson draw. The seed contains no sample-label information.
  for (int ieff = 0; ieff < kNEfficiency; ++ieff) {
    std::mt19937_64 random(splitmix64(
        0x4254474f63000000ULL ^
        (static_cast<std::uint64_t>(static_cast<std::uint32_t>(eventID)) << 8U) ^
        static_cast<std::uint64_t>(ieff)));
    int totalNpe = 0, maxNpe = 0;
    int nsegGe1 = 0, nsegGe2 = 0, nsegGe3 = 0;
    int matchedNpe = 0, matchedMaxNpe = 0;
    int matchedGe1 = 0, matchedGe2 = 0, matchedGe3 = 0;
    for (int i = 0; i < m.nSegTLC; ++i) {
      std::poisson_distribution<int> poisson(
          std::max(0.0, (*tlcExpected)[i] * kEfficiencies[ieff]));
      const int npe = poisson(random);
      totalNpe += npe;
      maxNpe = std::max(maxNpe, npe);
      if (npe >= 1) ++nsegGe1;
      if (npe >= 2) ++nsegGe2;
      if (npe >= 3) ++nsegGe3;
      if (matchedMask[i]) {
        matchedNpe += npe;
        matchedMaxNpe = std::max(matchedMaxNpe, npe);
        if (npe >= 1) ++matchedGe1;
        if (npe >= 2) ++matchedGe2;
        if (npe >= 3) ++matchedGe3;
      }
    }
    f[poissonFeature(ieff, 0)] = totalNpe;
    f[poissonFeature(ieff, 1)] = maxNpe;
    f[poissonFeature(ieff, 2)] = nsegGe1;
    f[poissonFeature(ieff, 3)] = nsegGe2;
    f[poissonFeature(ieff, 4)] = nsegGe3;
    f[poissonFeature(ieff, 5)] = matchedNpe;
    f[poissonFeature(ieff, 6)] =
        totalNpe > 0 ? static_cast<double>(matchedNpe) / totalNpe : 0.0;
    f[poissonFeature(ieff, 7)] = matchedMaxNpe;
    f[poissonFeature(ieff, 8)] = matchedGe1;
    f[poissonFeature(ieff, 9)] = matchedGe2;
    f[poissonFeature(ieff, 10)] = matchedGe3;
    f[poissonFeature(ieff, 11)] = nsegGe1 > 0 ? 1.0f : 0.0f;
    f[poissonFeature(ieff, 12)] = nsegGe2 > 0 ? 1.0f : 0.0f;
    f[poissonFeature(ieff, 13)] = nsegGe3 > 0 ? 1.0f : 0.0f;
    f[poissonFeature(ieff, 14)] = matchedGe1 > 0 ? 1.0f : 0.0f;
    f[poissonFeature(ieff, 15)] = matchedGe2 > 0 ? 1.0f : 0.0f;
    f[poissonFeature(ieff, 16)] = matchedGe3 > 0 ? 1.0f : 0.0f;
  }
}

void writeU32(std::ofstream &out, std::uint32_t v) {
  out.write(reinterpret_cast<const char *>(&v), sizeof(v));
}
void writeI32(std::ofstream &out, std::int32_t v) {
  out.write(reinterpret_cast<const char *>(&v), sizeof(v));
}
void writeF64(std::ofstream &out, double v) {
  out.write(reinterpret_cast<const char *>(&v), sizeof(v));
}

void writeManifest(const std::string &path, const std::string &input,
                   const Meta &m, double threshold, std::size_t nrow,
                   std::size_t missingTH, std::size_t missingTLC,
                   std::size_t extraTH, std::size_t extraTLC) {
  std::ofstream out(path + ".json");
  out << std::setprecision(15);
  out << "{\n"
      << "  \"format\": \"BGO2\",\n"
      << "  \"version\": " << kVersion << ",\n"
      << "  \"input\": \"" << input << "\",\n"
      << "  \"nrow\": " << nrow << ",\n"
      << "  \"ncol\": " << kNFeature << ",\n"
      << "  \"nLayer\": " << m.nLayer << ",\n"
      << "  \"nSector\": " << m.nSector << ",\n"
      << "  \"segmentationMode\": " << m.segmentationMode << ",\n"
      << "  \"geometryMode\": " << m.geometryMode << ",\n"
      << "  \"physicsFlag\": " << m.physicsFlag << ",\n"
      << "  \"thetaMin_deg\": " << m.thetaMinDeg << ",\n"
      << "  \"thetaMax_deg\": " << m.thetaMaxDeg << ",\n"
      << "  \"rMin_cm\": " << m.rMinCm << ",\n"
      << "  \"thickness_cm\": " << m.thicknessCm << ",\n"
      << "  \"bgoZOffset_cm\": " << m.bgoZOffsetCm << ",\n"
      << "  \"nSegTH\": " << m.nSegTH << ",\n"
      << "  \"nSegTLC\": " << m.nSegTLC << ",\n"
      << "  \"readoutMode\": \""
      << (m.finalReadoutMetadata ? "strict-timing" : "legacy-readout")
      << "\",\n"
      << "  \"finalReadoutMetadata\": "
      << (m.finalReadoutMetadata ? "true" : "false") << ",\n"
      << "  \"thRMin_mm\": " << m.thRMinMm << ",\n"
      << "  \"thRMax_mm\": " << m.thRMaxMm << ",\n"
      << "  \"thBarZMin_mm\": " << m.thBarZMinMm << ",\n"
      << "  \"thBarZMax_mm\": " << m.thBarZMaxMm << ",\n"
      << "  \"thEffectiveLightSpeed_mm_per_ns\": "
      << m.thEffectiveLightSpeedMmPerNs << ",\n"
      << "  \"thTimingSmearingAppliedInSource\": "
      << m.thTimingSmearingApplied << ",\n"
      << "  \"threshold_MeV\": " << threshold << ",\n"
      << "  \"missingTH\": " << missingTH << ",\n"
      << "  \"missingTLC\": " << missingTLC << ",\n"
      << "  \"extraTH\": " << extraTH << ",\n"
      << "  \"extraTLC\": " << extraTLC << ",\n"
      << "  \"eventJoin\": \"eventID\",\n"
      << "  \"poissonResponse\": {"
         "\"efficiencies\": [0.01, 0.02, 0.05, 0.10], "
         "\"npeThresholds\": [1, 2, 3], "
         "\"primaryHitDefinition\": \"any of 30 segments has Npe >= fixed threshold\", "
         "\"seed\": \"splitmix64(eventID, efficiencyIndex)\"},\n"
      << "  \"timingResponse\": ";
  if (m.finalReadoutMetadata) {
    out << "{\"perEndGaussianSigma_ns\": [0.0, 0.2, 0.5, 1.0], "
           "\"independentEnds\": true, \"vEffective_mm_per_ns\": "
        << m.thEffectiveLightSpeedMmPerNs
        << ", \"seed\": \"splitmix64(eventID, THsegment, sigmaIndex)\", "
           "\"status\": \"provisional_no_calibration\"},\n";
  } else {
    out << "null,\n";
  }
  out << "  \"classifierInputs\": [\"calarr.dE_MeV\", "
         "\"evt.EdepPC_MeV\", \"evt.EdepPCDown_MeV\", "
         "\"evt.EdepPCUp_MeV\", \"th.dE_MeV\", "
         "\"th.time_ns\", ";
  if (m.finalReadoutMetadata) out << "\"th.zReco_mm\", ";
  out << "\"tlc.cherenkovExpectedPhotons\", "
         "\"tlc.cherenkovTime_ns\"],\n"
      << "  \"forbiddenInputsRead\": [],\n"
      << "  \"features\": [\n";
  for (int i = 0; i < kNFeature; ++i)
    out << "    \"" << kFeatureNames[i] << "\"" << (i + 1 == kNFeature ? "\n" : ",\n");
  out << "  ]\n}\n";
}

bool extract(const std::string &input, const std::string &output,
             double threshold, bool allowMissing, bool allowNonFlag4,
             bool allowLegacyHodo) {
  TFile file(input.c_str(), "READ");
  if (file.IsZombie()) throw std::runtime_error("cannot open " + input);
  auto *cal = dynamic_cast<TTree *>(file.Get("calarr"));
  auto *evt = dynamic_cast<TTree *>(file.Get("evt"));
  auto *metaTree = dynamic_cast<TTree *>(file.Get("runmeta"));
  auto *th = dynamic_cast<TTree *>(file.Get("th"));
  auto *tlc = dynamic_cast<TTree *>(file.Get("tlc"));
  if (!cal || !evt || !metaTree || !th || !tlc)
    throw std::runtime_error("required tree missing (calarr/evt/runmeta/th/tlc)");
  const Meta meta = readMeta(metaTree, allowLegacyHodo);
  if (meta.physicsFlag != 4 && !allowNonFlag4)
    throw std::runtime_error("physicsFlag is not 4");
  const Geometry geometry(meta);

  const auto evtIndex = buildEventIndex(evt);
  const auto thIndex = buildEventIndex(th);
  const auto tlcIndex = buildEventIndex(tlc);

  requireBranch(cal, "eventID");
  requireBranch(cal, "dE_MeV");
  cal->SetBranchStatus("*", 0);
  cal->SetBranchStatus("eventID", 1);
  cal->SetBranchStatus("dE_MeV", 1);
  int calEventID = -1;
  std::vector<double> *calE = nullptr;
  cal->SetBranchAddress("eventID", &calEventID);
  cal->SetBranchAddress("dE_MeV", &calE);

  requireBranch(evt, "EdepPC_MeV");
  requireBranch(evt, "EdepPCDown_MeV");
  requireBranch(evt, "EdepPCUp_MeV");
  evt->SetBranchStatus("*", 0);
  evt->SetBranchStatus("eventID", 1);
  evt->SetBranchStatus("EdepPC_MeV", 1);
  evt->SetBranchStatus("EdepPCDown_MeV", 1);
  evt->SetBranchStatus("EdepPCUp_MeV", 1);
  int evtEventID = -1;
  double pcSumE = 0.0, pcDownE = 0.0, pcUpE = 0.0;
  evt->SetBranchAddress("eventID", &evtEventID);
  evt->SetBranchAddress("EdepPC_MeV", &pcSumE);
  evt->SetBranchAddress("EdepPCDown_MeV", &pcDownE);
  evt->SetBranchAddress("EdepPCUp_MeV", &pcUpE);

  requireBranch(th, "dE_MeV");
  requireBranch(th, "time_ns");
  const char *finalTHBranches[] = {"timeLeft_ns", "timeRight_ns",
                                   "timeLeftMinusRight_ns", "zReco_mm"};
  bool hasFinalTHBranches = true;
  for (const char *name : finalTHBranches)
    hasFinalTHBranches = hasFinalTHBranches && th->GetBranch(name);
  if (!hasFinalTHBranches && !allowLegacyHodo)
    throw std::runtime_error("final TH timing branches missing");
  const bool useTHZReco = meta.finalReadoutMetadata && hasFinalTHBranches;
  th->SetBranchStatus("*", 0);
  th->SetBranchStatus("eventID", 1);
  th->SetBranchStatus("dE_MeV", 1);
  th->SetBranchStatus("time_ns", 1);
  if (useTHZReco) th->SetBranchStatus("zReco_mm", 1);
  int thEventID = -1;
  std::vector<double> *thE = nullptr, *thTime = nullptr, *thZReco = nullptr;
  th->SetBranchAddress("eventID", &thEventID);
  th->SetBranchAddress("dE_MeV", &thE);
  th->SetBranchAddress("time_ns", &thTime);
  if (useTHZReco) th->SetBranchAddress("zReco_mm", &thZReco);

  requireBranch(tlc, "cherenkovExpectedPhotons");
  requireBranch(tlc, "cherenkovTime_ns");
  tlc->SetBranchStatus("*", 0);
  tlc->SetBranchStatus("eventID", 1);
  tlc->SetBranchStatus("cherenkovExpectedPhotons", 1);
  tlc->SetBranchStatus("cherenkovTime_ns", 1);
  int tlcEventID = -1;
  std::vector<double> *tlcExpected = nullptr, *tlcTime = nullptr;
  tlc->SetBranchAddress("eventID", &tlcEventID);
  tlc->SetBranchAddress("cherenkovExpectedPhotons", &tlcExpected);
  tlc->SetBranchAddress("cherenkovTime_ns", &tlcTime);

  std::unordered_set<int> calIDs;
  calIDs.reserve(static_cast<size_t>(cal->GetEntries() * 1.3 + 1));
  std::vector<float> all;
  all.reserve(static_cast<size_t>(cal->GetEntries()) * kNFeature);
  std::size_t missingTH = 0, missingTLC = 0;
  for (Long64_t i = 0; i < cal->GetEntries(); ++i) {
    cal->GetEntry(i);
    if (!calIDs.insert(calEventID).second)
      throw std::runtime_error("duplicate eventID in calarr: " +
                               std::to_string(calEventID));
    if (!calE || static_cast<int>(calE->size()) != geometry.NCell())
      throw std::runtime_error("calarr vector size mismatch for event " +
                               std::to_string(calEventID));
    auto f = makeBgoFeatures(geometry, *calE, calEventID, threshold);

    const auto eit = evtIndex.find(calEventID);
    if (eit == evtIndex.end())
      throw std::runtime_error("evt eventID is missing from calarr join: " +
                               std::to_string(calEventID));
    evt->GetEntry(eit->second);
    if (evtEventID != calEventID || !std::isfinite(pcSumE) ||
        !std::isfinite(pcDownE) || !std::isfinite(pcUpE) ||
        pcSumE < 0.0 || pcDownE < 0.0 || pcUpE < 0.0 ||
        std::abs(pcSumE - pcDownE - pcUpE) >
            1.0e-8 * std::max(1.0, pcSumE))
      throw std::runtime_error("evt PC join/value mismatch for event " +
                               std::to_string(calEventID));
    f[kPCSumE] = static_cast<float>(pcSumE);
    f[kPCDownE] = static_cast<float>(pcDownE);
    f[kPCUpE] = static_cast<float>(pcUpE);

    const std::vector<double> *thisTHE = nullptr, *thisTHTime = nullptr;
    const std::vector<double> *thisTHZReco = nullptr;
    const std::vector<double> *thisTLCExpected = nullptr, *thisTLCTime = nullptr;
    const auto hit = thIndex.find(calEventID);
    if (hit == thIndex.end()) {
      ++missingTH;
    } else {
      th->GetEntry(hit->second);
      if (thEventID != calEventID || !thE || !thTime ||
          static_cast<int>(thE->size()) != meta.nSegTH ||
          static_cast<int>(thTime->size()) != meta.nSegTH ||
          (useTHZReco &&
           (!thZReco || static_cast<int>(thZReco->size()) != meta.nSegTH)))
        throw std::runtime_error("TH join/vector mismatch for event " +
                                 std::to_string(calEventID));
      thisTHE = thE;
      thisTHTime = thTime;
      thisTHZReco = useTHZReco ? thZReco : nullptr;
    }
    const auto lit = tlcIndex.find(calEventID);
    if (lit == tlcIndex.end()) {
      ++missingTLC;
    } else {
      tlc->GetEntry(lit->second);
      if (tlcEventID != calEventID || !tlcExpected || !tlcTime ||
          static_cast<int>(tlcExpected->size()) != meta.nSegTLC ||
          static_cast<int>(tlcTime->size()) != meta.nSegTLC)
        throw std::runtime_error("TLC join/vector mismatch for event " +
                                 std::to_string(calEventID));
      thisTLCExpected = tlcExpected;
      thisTLCTime = tlcTime;
    }
    addHodoscopeFeatures(f, meta, calEventID, thisTHE, thisTHTime,
                         thisTHZReco,
                         thisTLCExpected, thisTLCTime);
    all.insert(all.end(), f.begin(), f.end());
  }
  std::size_t extraTH = 0, extraTLC = 0;
  for (const auto &item : evtIndex)
    if (calIDs.count(item.first) == 0)
      throw std::runtime_error("evt contains eventID absent from calarr: " +
                               std::to_string(item.first));
  for (const auto &item : thIndex)
    if (calIDs.count(item.first) == 0) ++extraTH;
  for (const auto &item : tlcIndex)
    if (calIDs.count(item.first) == 0) ++extraTLC;
  if (!allowMissing &&
      (missingTH != 0 || missingTLC != 0 || extraTH != 0 || extraTLC != 0))
    throw std::runtime_error(
        "TH/TLC eventID set differs from calarr; rerun with --allow-missing only if intentional");

  std::ofstream out(output, std::ios::binary);
  if (!out) throw std::runtime_error("cannot create " + output);
  writeU32(out, kMagic);
  writeU32(out, kVersion);
  writeU32(out, static_cast<std::uint32_t>(cal->GetEntries()));
  writeU32(out, kNFeature);
  writeI32(out, meta.nLayer);
  writeI32(out, meta.nSector);
  writeI32(out, meta.segmentationMode);
  writeI32(out, meta.physicsFlag);
  writeI32(out, meta.nSegTH);
  writeI32(out, meta.nSegTLC);
  writeF64(out, meta.thetaMinDeg);
  writeF64(out, meta.thetaMaxDeg);
  writeF64(out, threshold);
  out.write(reinterpret_cast<const char *>(all.data()), all.size() * sizeof(float));
  if (!out) throw std::runtime_error("failed while writing " + output);
  writeManifest(output, input, meta, threshold, cal->GetEntries(),
                missingTH, missingTLC, extraTH, extraTLC);
  std::cout << input << " -> " << output
            << " rows=" << cal->GetEntries() << " cols=" << kNFeature
            << " geometry=" << meta.nLayer << "x" << meta.nSector
            << " segmentationMode=" << meta.segmentationMode
            << " physicsFlag=" << meta.physicsFlag
            << " TH/TLC=" << meta.nSegTH << "/" << meta.nSegTLC
            << " finalReadout="
            << useTHZReco
            << " missing=" << missingTH << "/" << missingTLC
            << " extra=" << extraTH << "/" << extraTLC
            << " threshold=" << threshold << " MeV\n";
  return true;
}

void selfTestGeometry() {
  Meta equal;
  equal.nLayer = 10;
  equal.nSector = 20;
  equal.segmentationMode = 1;
  equal.thetaMinDeg = 5.666;
  equal.thetaMaxDeg = 170.302;
  equal.rMinCm = 30.0;
  equal.thicknessCm = 20.0;
  equal.nSegTH = 30;
  equal.nSegTLC = 30;
  const Geometry g(equal);
  if (static_cast<int>(g.directions.size()) != equal.nLayer * equal.nSector)
    throw std::runtime_error("equal-solid geometry size self-test failed");
  const double cmin = std::cos(equal.thetaMinDeg * kPi / 180.0);
  const double cmax = std::cos(equal.thetaMaxDeg * kPi / 180.0);
  for (int layer = 0; layer < equal.nLayer; ++layer) {
    const double expectedZ =
        cmin + (cmax - cmin) * (layer + 0.5) / equal.nLayer;
    for (int sector = 0; sector < equal.nSector; ++sector) {
      const auto &u = g.directions[g.Index(layer, sector)];
      const double norm =
          std::sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
      if (std::abs(norm - 1.0) > 1e-12 ||
          std::abs(u[2] - expectedZ) > 1e-12)
        throw std::runtime_error("equal-solid direction self-test failed");
    }
  }
  Meta uniform = equal;
  uniform.nLayer = 15;
  uniform.nSector = 15;
  uniform.segmentationMode = 0;
  const Geometry ug(uniform);
  const double expectedTheta =
      (uniform.thetaMinDeg +
       0.5 * (uniform.thetaMaxDeg - uniform.thetaMinDeg) / uniform.nLayer) *
      kPi / 180.0;
  if (std::abs(ug.directions[0][2] - std::cos(expectedTheta)) > 1e-12)
    throw std::runtime_error("uniform-theta direction self-test failed");

  Meta egg;
  egg.nLayer = 31;
  egg.nSector = 60;
  egg.segmentationMode = 2;
  egg.geometryMode = 2;
  egg.thetaMinDeg = BGOeggGeometry::BuildRings(31).front().thetaLowDeg;
  egg.thetaMaxDeg = BGOeggGeometry::BuildRings(31).back().thetaHighDeg;
  egg.rMinCm = 20.0;
  egg.thicknessCm = 22.0;
  egg.nSegTH = 30;
  egg.nSegTLC = 30;
  const Geometry eg(egg);
  if (static_cast<int>(eg.directions.size()) != 1860)
    throw std::runtime_error("BGOegg published geometry size self-test failed");
  if (!(std::acos(eg.directions.front()[2]) * 180.0 / kPi > 5.0 &&
        std::acos(eg.directions.front()[2]) * 180.0 / kPi < 9.0))
    throw std::runtime_error("BGOegg first centroid direction self-test failed");
  Meta shifted = egg;
  shifted.bgoZOffsetCm = -10.0;
  const Geometry shiftedGeometry(shifted);
  if (!(shiftedGeometry.directions.front()[2] < eg.directions.front()[2]))
    throw std::runtime_error("BGOegg z-offset direction self-test failed");
  std::cout << "geometry self-test ok: equal-solid, uniform, and BGOegg published 31-layer\n";
}

}  // namespace

int main(int argc, char **argv) {
  if (argc == 2 && std::string(argv[1]) == "--self-test-geometry") {
    try {
      selfTestGeometry();
      return 0;
    } catch (const std::exception &e) {
      std::cerr << "ERROR: " << e.what() << '\n';
      return 1;
    }
  }
  if (argc < 4) {
    std::cerr << "usage: bgo_extract_features_v2 input.root output.bgo2 threshold_MeV"
                 " [--allow-missing] [--allow-nonflag4] [--allow-legacy-hodo]\n";
    return 2;
  }
  bool allowMissing = false, allowNonFlag4 = false, allowLegacyHodo = false;
  for (int i = 4; i < argc; ++i) {
    const std::string arg(argv[i]);
    if (arg == "--allow-missing") allowMissing = true;
    else if (arg == "--allow-nonflag4") allowNonFlag4 = true;
    else if (arg == "--allow-legacy-hodo") allowLegacyHodo = true;
    else {
      std::cerr << "unknown option: " << arg << '\n';
      return 2;
    }
  }
  try {
    return extract(argv[1], argv[2], std::stod(argv[3]),
                   allowMissing, allowNonFlag4, allowLegacyHodo) ? 0 : 1;
  } catch (const std::exception &e) {
    std::cerr << "ERROR: " << e.what() << '\n';
    return 1;
  }
}
