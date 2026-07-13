#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

namespace {

void requireBranch(TTree *tree, const char *name) {
  if (!tree || !tree->GetBranch(name))
    throw std::runtime_error(std::string("missing branch ") +
                             (tree ? tree->GetName() : "<null>") + "." + name);
}

std::unordered_set<int> eventIDs(TTree *tree) {
  requireBranch(tree, "eventID");
  tree->SetBranchStatus("*", 0);
  tree->SetBranchStatus("eventID", 1);
  int eventID = -1;
  tree->SetBranchAddress("eventID", &eventID);
  std::unordered_set<int> result;
  for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry) {
    tree->GetEntry(entry);
    if (!result.insert(eventID).second)
      throw std::runtime_error(
          std::string("duplicate eventID in ") + tree->GetName());
  }
  return result;
}

void requireFinite(const std::vector<double> &values,
                   const char *tree, const char *branch, int eventID) {
  for (double value : values)
    if (!std::isfinite(value))
      throw std::runtime_error(
          std::string("non-finite ") + tree + "." + branch +
          " eventID=" + std::to_string(eventID));
}

void validate(const std::string &path) {
  TFile file(path.c_str(), "READ");
  if (file.IsZombie()) throw std::runtime_error("cannot open " + path);
  auto *meta = dynamic_cast<TTree *>(file.Get("runmeta"));
  auto *cal = dynamic_cast<TTree *>(file.Get("calarr"));
  auto *th = dynamic_cast<TTree *>(file.Get("th"));
  auto *tlc = dynamic_cast<TTree *>(file.Get("tlc"));
  if (!meta || !cal || !th || !tlc)
    throw std::runtime_error("required tree missing");
  if (meta->GetEntries() != 1)
    throw std::runtime_error("runmeta entries != 1");

  int nLayer = 0, nSector = 0, segmentationMode = -1;
  int physicsFlag = -1, nSegTH = 0, nSegTLC = 0;
  double thRMin = 0.0, thRMax = 0.0, thZMin = 0.0, thZMax = 0.0;
  double thLightSpeed = 0.0;
  for (const char *name : {"nLayer", "nSector", "segmentationMode",
                           "physicsFlag", "nSegTH", "nSegTLC",
                           "thRMin_mm", "thRMax_mm", "thBarZMin_mm",
                           "thBarZMax_mm",
                           "thEffectiveLightSpeed_mm_per_ns",
                           "thTimingSmearingApplied", "thTimingModel"})
    requireBranch(meta, name);
  meta->SetBranchStatus("*", 0);
  for (const char *name : {"nLayer", "nSector", "segmentationMode",
                           "physicsFlag", "nSegTH", "nSegTLC",
                           "thRMin_mm", "thRMax_mm", "thBarZMin_mm",
                           "thBarZMax_mm",
                           "thEffectiveLightSpeed_mm_per_ns"})
    meta->SetBranchStatus(name, 1);
  meta->SetBranchAddress("nLayer", &nLayer);
  meta->SetBranchAddress("nSector", &nSector);
  meta->SetBranchAddress("segmentationMode", &segmentationMode);
  meta->SetBranchAddress("physicsFlag", &physicsFlag);
  meta->SetBranchAddress("nSegTH", &nSegTH);
  meta->SetBranchAddress("nSegTLC", &nSegTLC);
  meta->SetBranchAddress("thRMin_mm", &thRMin);
  meta->SetBranchAddress("thRMax_mm", &thRMax);
  meta->SetBranchAddress("thBarZMin_mm", &thZMin);
  meta->SetBranchAddress("thBarZMax_mm", &thZMax);
  meta->SetBranchAddress("thEffectiveLightSpeed_mm_per_ns", &thLightSpeed);
  meta->GetEntry(0);
  if (nLayer <= 0 || nSector <= 0 || nSegTH <= 0 || nSegTLC <= 0 ||
      (segmentationMode != 0 && segmentationMode != 1))
    throw std::runtime_error("invalid runmeta values");
  if (!(thRMin > 0.0 && thRMin < thRMax && thZMin < thZMax &&
        thLightSpeed > 0.0))
    throw std::runtime_error("invalid TH runmeta geometry/timing values");
  if (physicsFlag != 4)
    throw std::runtime_error("physicsFlag != 4");

  const auto calIDs = eventIDs(cal);
  const auto thIDs = eventIDs(th);
  const auto tlcIDs = eventIDs(tlc);
  if (calIDs != thIDs || calIDs != tlcIDs)
    throw std::runtime_error("calarr/TH/TLC eventID sets differ");

  requireBranch(cal, "dE_MeV");
  cal->SetBranchStatus("*", 0);
  cal->SetBranchStatus("eventID", 1);
  cal->SetBranchStatus("dE_MeV", 1);
  int eventID = -1;
  std::vector<double> *calE = nullptr;
  cal->SetBranchAddress("eventID", &eventID);
  cal->SetBranchAddress("dE_MeV", &calE);
  for (Long64_t entry = 0; entry < cal->GetEntries(); ++entry) {
    cal->GetEntry(entry);
    if (!calE || static_cast<int>(calE->size()) != nLayer * nSector)
      throw std::runtime_error("calarr vector length mismatch");
    requireFinite(*calE, "calarr", "dE_MeV", eventID);
  }

  const char *thBranches[] = {
      "dE_MeV", "time_ns", "timeLeft_ns", "timeRight_ns",
      "timeLeftMinusRight_ns", "zReco_mm", "chargedPath_truth_mm"};
  for (const char *name : thBranches)
    requireBranch(th, name);
  th->SetBranchStatus("*", 0);
  th->SetBranchStatus("eventID", 1);
  for (const char *name : thBranches) th->SetBranchStatus(name, 1);
  std::vector<double> *thE = nullptr, *thTime = nullptr, *thPath = nullptr;
  std::vector<double> *thLeft = nullptr, *thRight = nullptr;
  std::vector<double> *thDt = nullptr, *thZReco = nullptr;
  th->SetBranchAddress("eventID", &eventID);
  th->SetBranchAddress("dE_MeV", &thE);
  th->SetBranchAddress("time_ns", &thTime);
  th->SetBranchAddress("timeLeft_ns", &thLeft);
  th->SetBranchAddress("timeRight_ns", &thRight);
  th->SetBranchAddress("timeLeftMinusRight_ns", &thDt);
  th->SetBranchAddress("zReco_mm", &thZReco);
  th->SetBranchAddress("chargedPath_truth_mm", &thPath);
  long long validTHZSegments = 0;
  double maxTHZFormulaResidual = 0.0;
  for (Long64_t entry = 0; entry < th->GetEntries(); ++entry) {
    th->GetEntry(entry);
    if (!thE || !thTime || !thLeft || !thRight || !thDt || !thZReco ||
        !thPath ||
        static_cast<int>(thE->size()) != nSegTH ||
        static_cast<int>(thTime->size()) != nSegTH ||
        static_cast<int>(thLeft->size()) != nSegTH ||
        static_cast<int>(thRight->size()) != nSegTH ||
        static_cast<int>(thDt->size()) != nSegTH ||
        static_cast<int>(thZReco->size()) != nSegTH ||
        static_cast<int>(thPath->size()) != nSegTH)
      throw std::runtime_error("TH vector length mismatch");
    requireFinite(*thE, "th", "dE_MeV", eventID);
    requireFinite(*thTime, "th", "time_ns", eventID);
    requireFinite(*thLeft, "th", "timeLeft_ns", eventID);
    requireFinite(*thRight, "th", "timeRight_ns", eventID);
    requireFinite(*thDt, "th", "timeLeftMinusRight_ns", eventID);
    requireFinite(*thZReco, "th", "zReco_mm", eventID);
    requireFinite(*thPath, "th", "chargedPath_truth_mm", eventID);
    for (int segment = 0; segment < nSegTH; ++segment) {
      if ((*thPath)[segment] < 0.0)
        throw std::runtime_error("negative TH truth path");
      const bool leftValid = (*thLeft)[segment] >= 0.0;
      const bool rightValid = (*thRight)[segment] >= 0.0;
      if (leftValid != rightValid)
        throw std::runtime_error("one-sided TH timing validity");
      if (!leftValid) {
        if ((*thDt)[segment] != -9999.0 || (*thZReco)[segment] != -9999.0)
          throw std::runtime_error("TH timing sentinel mismatch");
        continue;
      }
      ++validTHZSegments;
      const double dtExpected = (*thLeft)[segment] - (*thRight)[segment];
      const double zExpected = 0.5 * thLightSpeed * dtExpected;
      maxTHZFormulaResidual = std::max(
          maxTHZFormulaResidual,
          std::max(std::abs((*thDt)[segment] - dtExpected),
                   std::abs((*thZReco)[segment] - zExpected)));
      if (std::abs((*thDt)[segment] - dtExpected) > 1e-9 ||
          std::abs((*thZReco)[segment] - zExpected) > 1e-8 ||
          (*thZReco)[segment] < thZMin - 1e-8 ||
          (*thZReco)[segment] > thZMax + 1e-8)
        throw std::runtime_error("TH reconstructed-z consistency failure");
    }
  }

  const char *tlcBranches[] = {
      "dE_truth_MeV", "cherenkovTime_ns", "chargedPath_truth_mm",
      "cherenkovPath_mm", "cherenkovExpectedPhotons"};
  for (const char *name : tlcBranches) requireBranch(tlc, name);
  tlc->SetBranchStatus("*", 0);
  tlc->SetBranchStatus("eventID", 1);
  for (const char *name : tlcBranches) tlc->SetBranchStatus(name, 1);
  std::vector<double> *tlcE = nullptr, *tlcTime = nullptr;
  std::vector<double> *truthPath = nullptr, *cherenkovPath = nullptr;
  std::vector<double> *expectedPhotons = nullptr;
  tlc->SetBranchAddress("eventID", &eventID);
  tlc->SetBranchAddress("dE_truth_MeV", &tlcE);
  tlc->SetBranchAddress("cherenkovTime_ns", &tlcTime);
  tlc->SetBranchAddress("chargedPath_truth_mm", &truthPath);
  tlc->SetBranchAddress("cherenkovPath_mm", &cherenkovPath);
  tlc->SetBranchAddress("cherenkovExpectedPhotons", &expectedPhotons);

  double minYield = 1e300, maxYield = 0.0, maxPathExcess = -1e300;
  long long positiveYieldSegments = 0;
  for (Long64_t entry = 0; entry < tlc->GetEntries(); ++entry) {
    tlc->GetEntry(entry);
    for (auto *values : {tlcE, tlcTime, truthPath,
                         cherenkovPath, expectedPhotons})
      if (!values || static_cast<int>(values->size()) != nSegTLC)
        throw std::runtime_error("TLC vector length mismatch");
    requireFinite(*tlcE, "tlc", "dE_truth_MeV", eventID);
    requireFinite(*tlcTime, "tlc", "cherenkovTime_ns", eventID);
    requireFinite(*truthPath, "tlc", "chargedPath_truth_mm", eventID);
    requireFinite(*cherenkovPath, "tlc", "cherenkovPath_mm", eventID);
    requireFinite(*expectedPhotons, "tlc", "cherenkovExpectedPhotons", eventID);
    for (int segment = 0; segment < nSegTLC; ++segment) {
      if ((*truthPath)[segment] < 0.0 || (*cherenkovPath)[segment] < 0.0)
        throw std::runtime_error("negative TLC path");
      const double excess =
          (*cherenkovPath)[segment] - (*truthPath)[segment];
      maxPathExcess = std::max(maxPathExcess, excess);
      if (excess > 1e-9)
        throw std::runtime_error(
            "cherenkovPath_mm > chargedPath_truth_mm eventID=" +
            std::to_string(eventID));
      minYield = std::min(minYield, (*expectedPhotons)[segment]);
      maxYield = std::max(maxYield, (*expectedPhotons)[segment]);
      if ((*expectedPhotons)[segment] < 0.0)
        throw std::runtime_error("negative Cherenkov expected yield");
      if ((*expectedPhotons)[segment] > 0.0) ++positiveYieldSegments;
    }
  }
  std::cout << "VALID " << path
            << " geometry=" << nLayer << "x" << nSector
            << " segmentationMode=" << segmentationMode
            << " events=" << calIDs.size()
            << " vectors=" << nLayer * nSector << "/" << nSegTH
            << "/" << nSegTLC
            << " eventID_join=strict"
            << " TH[rMin,rMax]=" << thRMin << "," << thRMax
            << " validTHZSegments=" << validTHZSegments
            << " maxTHZFormulaResidual=" << maxTHZFormulaResidual
            << " max(cherenkovPath-truthPath)=" << maxPathExcess
            << " expectedYield[min,max]=" << minYield << "," << maxYield
            << " positiveYieldSegments=" << positiveYieldSegments
            << " [QA-only truth path read; classifier does not read truth]\n";
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    std::cerr << "usage: validate_hodo_root_schema file.root [file.root ...]\n";
    return 2;
  }
  try {
    for (int i = 1; i < argc; ++i) validate(argv[i]);
    return 0;
  } catch (const std::exception &error) {
    std::cerr << "ERROR: " << error.what() << '\n';
    return 1;
  }
}