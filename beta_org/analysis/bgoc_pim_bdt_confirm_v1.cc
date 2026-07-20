#include <TMVA/Reader.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

constexpr std::uint32_t kMagic = 0x42474f32U;
constexpr std::uint32_t kVersion = 2U;
constexpr double kThreshold = 0.480194806111;
constexpr const char *kMethod = "BDTG_D3_M1";

struct Bgo2Header {
  std::uint32_t magic;
  std::uint32_t version;
  std::uint32_t nrow;
  std::uint32_t ncol;
  std::int32_t nLayer;
  std::int32_t nSector;
  std::int32_t segmentationMode;
  std::int32_t physicsFlag;
  std::int32_t nSegTH;
  std::int32_t nSegTLC;
  double thetaMinDeg;
  double thetaMaxDeg;
  double thresholdMeV;
};

static_assert(sizeof(Bgo2Header) == 64, "unexpected BGO2 header layout");

const std::array<std::string, 23> kFeatures = {
    "maxE",          "sumE",          "nHit",          "nCl4",
    "nCl8",          "cl1Sum4",       "cl2Sum4",       "cl1Size4",
    "cl1MaxFrac4",   "cl1RmsDeg4",    "allRmsDeg",     "cl2OverCl1",
    "isolatedEFrac4", "seed15NCl4",    "seed15Cl1Sum4", "seed15Cl1Size4",
    "e2OverE1",      "nHit2T",        "sum2T",         "local3x3Frac",
    "meanHitE",      "thetaSpan",     "phiSpan"};

struct Interval {
  double low{};
  double high{};
};

struct Metric {
  std::uint64_t numerator{};
  std::uint64_t denominator{};
  double value{};
  Interval interval{};
};

Interval wilson(std::uint64_t numerator, std::uint64_t denominator) {
  constexpr double z = 1.959963984540054;
  const double n = static_cast<double>(denominator);
  const double p = static_cast<double>(numerator) / n;
  const double denominatorTerm = 1.0 + z * z / n;
  const double centre = (p + z * z / (2.0 * n)) / denominatorTerm;
  const double halfWidth =
      z * std::sqrt((p * (1.0 - p) + z * z / (4.0 * n)) / n) /
      denominatorTerm;
  return {centre - halfWidth, centre + halfWidth};
}

Metric metric(std::uint64_t numerator, std::uint64_t denominator) {
  return {numerator, denominator,
          static_cast<double>(numerator) / static_cast<double>(denominator),
          wilson(numerator, denominator)};
}

template <class Callback>
void forEachBgo2Row(const fs::path &path, Callback callback) {
  std::ifstream input(path, std::ios::binary);
  if (!input) throw std::runtime_error("cannot open " + path.string());
  Bgo2Header header{};
  input.read(reinterpret_cast<char *>(&header), sizeof(header));
  if (!input || header.magic != kMagic || header.version != kVersion ||
      header.ncol < 24 || header.nLayer != 15 || header.nSector != 15 ||
      header.segmentationMode != 0 || header.physicsFlag != 4 ||
      std::abs(header.thresholdMeV - 3.0) > 1.0e-12) {
    throw std::runtime_error("incompatible BGO2 header in " + path.string());
  }
  std::vector<float> row(header.ncol);
  for (std::uint32_t entry = 0; entry < header.nrow; ++entry) {
    input.read(reinterpret_cast<char *>(row.data()),
               static_cast<std::streamsize>(sizeof(float) * row.size()));
    if (!input) throw std::runtime_error("truncated BGO2 row in " + path.string());
    callback(row);
  }
  if (input.peek() != std::ifstream::traits_type::eof()) {
    throw std::runtime_error("extra bytes in " + path.string());
  }
}

std::vector<double> evaluate(const fs::path &path, TMVA::Reader &reader,
                             std::array<float, 23> &values) {
  std::vector<double> scores;
  forEachBgo2Row(path, [&](const std::vector<float> &row) {
    for (std::size_t i = 0; i < kFeatures.size(); ++i) values[i] = row[i + 1];
    scores.push_back(reader.EvaluateMVA(kMethod));
  });
  return scores;
}

double auc(const std::vector<double> &electron,
           const std::vector<double> &pim) {
  struct Item {
    double score;
    bool isPim;
  };
  std::vector<Item> pooled;
  pooled.reserve(electron.size() + pim.size());
  for (double value : electron) pooled.push_back({value, false});
  for (double value : pim) pooled.push_back({value, true});
  std::sort(pooled.begin(), pooled.end(),
            [](const Item &a, const Item &b) { return a.score < b.score; });
  double pimRankSum = 0.0;
  std::size_t begin = 0;
  while (begin < pooled.size()) {
    std::size_t end = begin + 1;
    while (end < pooled.size() && pooled[end].score == pooled[begin].score) ++end;
    const double averageRank =
        0.5 * (static_cast<double>(begin + 1) + static_cast<double>(end));
    for (std::size_t i = begin; i < end; ++i) {
      if (pooled[i].isPim) pimRankSum += averageRank;
    }
    begin = end;
  }
  const double np = static_cast<double>(pim.size());
  const double ne = static_cast<double>(electron.size());
  const double uPim = pimRankSum - np * (np + 1.0) / 2.0;
  return 1.0 - uPim / (ne * np);
}

void writeMetric(std::ostream &out, const Metric &value) {
  out << "{\"numerator\":" << value.numerator
      << ",\"denominator\":" << value.denominator
      << ",\"value\":" << value.value << ",\"wilson_95ci\":["
      << value.interval.low << ',' << value.interval.high << "]}";
}

}  // namespace

int main(int argc, char **argv) {
  try {
    if (argc != 3) {
      std::cerr << "usage: bgoc_pim_bdt_confirm_v1 PROJECT INPUT_BGO2_DIR\n";
      return 2;
    }
    const fs::path project = fs::path(argv[1]);
    const fs::path inputDirectory = fs::path(argv[2]);
    const fs::path model =
        project / "analysis/models/bgoc_pim_bdtg_d3_m1_seed6302026.weights.xml";
    const fs::path result = project /
        "analysis/results/bgoegg_v1/bgoc_pim_bdt_confirm_s10302026.json";
    const fs::path rootOutput =
        project / "tmp/bgoc_pim_bdt_confirm_s10302026/scores.root";
    fs::create_directories(result.parent_path());
    fs::create_directories(rootOutput.parent_path());

    TMVA::Reader reader("!Color:!Silent");
    std::array<float, 23> values{};
    for (std::size_t i = 0; i < kFeatures.size(); ++i) {
      reader.AddVariable(kFeatures[i].c_str(), &values[i]);
    }
    reader.BookMVA(kMethod, model.c_str());

    const auto electron = evaluate(inputDirectory / "e.bgo2", reader, values);
    const auto pim = evaluate(inputDirectory / "pim.bgo2", reader, values);
    const auto pi0 = evaluate(inputDirectory / "pi0.bgo2", reader, values);
    if (electron.size() != 100000 || pim.size() != 100000 ||
        pi0.size() != 100000) {
      throw std::runtime_error("confirmation requires exactly 100000 events/species");
    }

    const auto count = [](const std::vector<double> &scores, auto predicate) {
      return static_cast<std::uint64_t>(
          std::count_if(scores.begin(), scores.end(), predicate));
    };
    const Metric electronKeep = metric(
        count(electron, [](double score) { return score > kThreshold; }),
        electron.size());
    const Metric pimReject = metric(
        count(pim, [](double score) { return score <= kThreshold; }), pim.size());
    const Metric pi0Reject = metric(
        count(pi0, [](double score) { return score <= kThreshold; }), pi0.size());

    std::ofstream json(result);
    if (!json) throw std::runtime_error("cannot create " + result.string());
    json << std::setprecision(12);
    json << "{\n  \"status\":\"frozen_independent_confirmation\",\n"
         << "  \"analysis\":\"bgoc_pim_bdt_confirm_v1\",\n"
         << "  \"objective\":\"pim reject > 90% and electron keep > 60%\",\n"
         << "  \"sample\":\"seed10302026, 100000 events/species, untouched before model and cut freeze\",\n"
         << "  \"model\":\"analysis/models/bgoc_pim_bdtg_d3_m1_seed6302026.weights.xml\",\n"
         << "  \"method\":\"" << kMethod << "\",\n"
         << "  \"threshold\":" << kThreshold << ",\n"
         << "  \"cut_rule\":\"accept electron if score > threshold\",\n"
         << "  \"pi0_role\":\"guard only; not used in training or cut selection\",\n"
         << "  \"features\":[";
    for (std::size_t i = 0; i < kFeatures.size(); ++i) {
      if (i) json << ',';
      json << '\"' << kFeatures[i] << '\"';
    }
    json << "],\n  \"metrics\":{\n    \"electron_keep\":";
    writeMetric(json, electronKeep);
    json << ",\n    \"pim_reject\":";
    writeMetric(json, pimReject);
    json << ",\n    \"pi0_reject\":";
    writeMetric(json, pi0Reject);
    json << ",\n    \"auc_e_pim\":" << auc(electron, pim)
         << "\n  },\n  \"passes_strict_target\":"
         << ((electronKeep.value > 0.60 && pimReject.value > 0.90) ? "true" : "false")
         << "\n}\n";

    TFile output(rootOutput.c_str(), "RECREATE");
    if (output.IsZombie()) {
      throw std::runtime_error("cannot create " + rootOutput.string());
    }
    TH1D hElectron("hScore_e", "electron;BDTG score;events", 160, -0.2, 0.8);
    TH1D hPim("hScore_pim", "pi-minus;BDTG score;events", 160, -0.2, 0.8);
    TH1D hPi0("hScore_pi0", "pi-zero;BDTG score;events", 160, -0.2, 0.8);
    for (double score : electron) hElectron.Fill(score);
    for (double score : pim) hPim.Fill(score);
    for (double score : pi0) hPi0.Fill(score);
    hElectron.Write();
    hPim.Write();
    hPi0.Write();

    constexpr int nRoc = 1001;
    TGraph roc(nRoc);
    roc.SetName("gRoc_e_pim");
    roc.SetTitle("electron efficiency vs pi-minus rejection;pi-minus rejection;electron efficiency");
    std::vector<double> ordered = pim;
    std::sort(ordered.begin(), ordered.end());
    for (int i = 0; i < nRoc; ++i) {
      const double reject = static_cast<double>(i) / (nRoc - 1);
      const std::size_t index = std::min<std::size_t>(
          ordered.size() - 1,
          static_cast<std::size_t>(std::floor(reject * (ordered.size() - 1))));
      const double threshold = ordered[index];
      const double keep = static_cast<double>(
          count(electron, [&](double score) { return score > threshold; })) /
          electron.size();
      roc.SetPoint(i, reject, keep);
    }
    roc.Write();
    output.Close();

    std::cout << std::fixed << std::setprecision(3)
              << "electron keep=" << 100.0 * electronKeep.value
              << "% pim reject=" << 100.0 * pimReject.value
              << "% pi0 reject=" << 100.0 * pi0Reject.value
              << "% pass="
              << (electronKeep.value > 0.60 && pimReject.value > 0.90 ? "yes" : "no")
              << '\n'
              << result << '\n' << rootOutput << '\n';
    return 0;
  } catch (const std::exception &error) {
    std::cerr << "ERROR: " << error.what() << '\n';
    return 1;
  }
}
