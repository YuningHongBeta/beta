#include <TMVA/DataLoader.h>
#include <TMVA/Factory.h>
#include <TMVA/Reader.h>
#include <TMVA/Tools.h>
#include <TMVA/Types.h>
#include <TFile.h>
#include <TRandom.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

constexpr std::uint32_t kMagic = 0x42474f32U;
constexpr std::uint32_t kVersion = 2U;
constexpr std::size_t kBgoFeatureCount = 23;
constexpr std::size_t kMaximumFeatureCount = 24;
constexpr std::size_t kPcUpColumn = 64;
constexpr const char *kMethod = "BDTG_D3_M1";
constexpr const char *kOptions =
    "!H:!V:NTrees=800:MinNodeSize=1.0%:BoostType=Grad:Shrinkage=0.03:"
    "UseBaggedBoost:BaggedSampleFraction=0.6:nCuts=40:MaxDepth=3";
constexpr double kPi0ToBetaFactor = 0.2 * 0.1 * 0.1 * 0.8 / 0.00048;
constexpr double kPimToBetaFactor = 0.4 * 0.001 * 0.045 / 0.00048;
constexpr double kTotalRequirement = 0.04;

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

const std::array<std::string, kMaximumFeatureCount> kFeatures = {
    "maxE",           "sumE",           "nHit",          "nCl4",
    "nCl8",           "cl1Sum4",        "cl2Sum4",       "cl1Size4",
    "cl1MaxFrac4",    "cl1RmsDeg4",     "allRmsDeg",     "cl2OverCl1",
    "isolatedEFrac4", "seed15NCl4",     "seed15Cl1Sum4", "seed15Cl1Size4",
    "e2OverE1",       "nHit2T",         "sum2T",         "local3x3Frac",
    "meanHitE",       "thetaSpan",      "phiSpan",       "pcUpE_MeV"};

struct GeometrySpec {
  std::string name;
  std::string prefix;
  int nLayer;
  int nSector;
  int segmentation;
};

const std::array<GeometrySpec, 2> kGeometries = {{
    {"BGOC", "bgoc", 15, 15, 0},
    {"BGOegg31", "bgoegg31", 31, 60, 2},
}};

struct OwnedTree {
  std::unique_ptr<TTree> tree;
  std::array<float, kMaximumFeatureCount> values{};

  OwnedTree(const std::string &name, std::size_t nFeature)
      : tree(std::make_unique<TTree>(name.c_str(), name.c_str())) {
    for (std::size_t index = 0; index < nFeature; ++index) {
      tree->Branch(kFeatures[index].c_str(), &values[index],
                   (kFeatures[index] + "/F").c_str());
    }
  }
};

struct SampleScores {
  std::vector<double> electron;
  std::vector<double> pim;
  std::vector<double> pi0;
};

struct Counts {
  std::uint64_t electronKeep{};
  std::uint64_t pimReject{};
  std::uint64_t pi0Reject{};
  std::uint64_t denominator{};
};

template <class Callback>
void forEachRow(const fs::path &path, const GeometrySpec &geometry,
                Callback callback) {
  std::ifstream input(path, std::ios::binary);
  if (!input) throw std::runtime_error("cannot open " + path.string());
  Bgo2Header header{};
  input.read(reinterpret_cast<char *>(&header), sizeof(header));
  if (!input || header.magic != kMagic || header.version != kVersion ||
      header.ncol <= kPcUpColumn || header.nLayer != geometry.nLayer ||
      header.nSector != geometry.nSector ||
      header.segmentationMode != geometry.segmentation ||
      header.physicsFlag != 4 || std::abs(header.thresholdMeV - 3.0) > 1e-12) {
    throw std::runtime_error("incompatible BGO2 header in " + path.string());
  }
  std::vector<float> row(header.ncol);
  for (std::uint32_t entry = 0; entry < header.nrow; ++entry) {
    input.read(reinterpret_cast<char *>(row.data()),
               static_cast<std::streamsize>(row.size() * sizeof(float)));
    if (!input) throw std::runtime_error("truncated row in " + path.string());
    const auto eventId = static_cast<std::int64_t>(std::llround(row[0]));
    if (static_cast<float>(eventId) != row[0]) {
      throw std::runtime_error("non-integral eventID in " + path.string());
    }
    callback(eventId, row);
  }
  if (input.peek() != std::ifstream::traits_type::eof()) {
    throw std::runtime_error("extra bytes in " + path.string());
  }
}

void copyFeatures(const std::vector<float> &row,
                  std::array<float, kMaximumFeatureCount> &values,
                  bool includePc) {
  for (std::size_t index = 0; index < kBgoFeatureCount; ++index) {
    values[index] = row[index + 1];
  }
  if (includePc) values[kBgoFeatureCount] = row[kPcUpColumn];
}

void fillTree(const fs::path &path, const GeometrySpec &geometry,
              bool includePc, OwnedTree &output) {
  forEachRow(path, geometry,
             [&](std::int64_t eventId, const std::vector<float> &row) {
    if ((eventId & 3) >= 2) return;
    copyFeatures(row, output.values, includePc);
    output.tree->Fill();
  });
}

std::vector<double> evaluateFile(
    const fs::path &path, const GeometrySpec &geometry, bool includePc,
    int split, TMVA::Reader &reader,
    std::array<float, kMaximumFeatureCount> &values) {
  std::vector<double> scores;
  forEachRow(path, geometry,
             [&](std::int64_t eventId, const std::vector<float> &row) {
    if (split >= 0 && (eventId & 3) != split) return;
    copyFeatures(row, values, includePc);
    scores.push_back(reader.EvaluateMVA(kMethod));
  });
  return scores;
}

SampleScores evaluateSample(
    const fs::path &directory, const GeometrySpec &geometry, bool includePc,
    int split, TMVA::Reader &reader,
    std::array<float, kMaximumFeatureCount> &values) {
  const auto path = [&](const std::string &species) {
    return directory /
        (geometry.prefix + "_balanced_pc3_h60_s7232026_" + species + ".bgo2");
  };
  return {evaluateFile(path("e"), geometry, includePc, split, reader, values),
          evaluateFile(path("pim"), geometry, includePc, split, reader, values),
          evaluateFile(path("pi0"), geometry, includePc, split, reader, values)};
}

double auc(const std::vector<double> &electron,
           const std::vector<double> &background) {
  struct Item { double score; bool background; };
  std::vector<Item> pooled;
  pooled.reserve(electron.size() + background.size());
  for (double score : electron) pooled.push_back({score, false});
  for (double score : background) pooled.push_back({score, true});
  std::sort(pooled.begin(), pooled.end(),
            [](const Item &left, const Item &right) {
              return left.score < right.score;
            });
  double backgroundRankSum = 0.0;
  std::size_t begin = 0;
  while (begin < pooled.size()) {
    std::size_t end = begin + 1;
    while (end < pooled.size() && pooled[end].score == pooled[begin].score) ++end;
    const double averageRank =
        0.5 * (static_cast<double>(begin + 1) + static_cast<double>(end));
    for (std::size_t index = begin; index < end; ++index) {
      if (pooled[index].background) backgroundRankSum += averageRank;
    }
    begin = end;
  }
  const double nb = static_cast<double>(background.size());
  const double ne = static_cast<double>(electron.size());
  const double uBackground = backgroundRankSum - nb * (nb + 1.0) / 2.0;
  return 1.0 - uBackground / (ne * nb);
}

double thresholdForElectronKeep(const std::vector<double> &scores,
                                double targetKeep) {
  std::vector<double> ordered = scores;
  std::sort(ordered.begin(), ordered.end());
  const std::size_t rejectCount = static_cast<std::size_t>(
      std::floor((1.0 - targetKeep) * ordered.size()));
  if (rejectCount == 0) return std::nextafter(ordered.front(), -INFINITY);
  if (rejectCount >= ordered.size()) return ordered.back();
  return 0.5 * (ordered[rejectCount - 1] + ordered[rejectCount]);
}

Counts scoreCounts(const SampleScores &scores, double threshold) {
  const auto count = [](const std::vector<double> &values, auto predicate) {
    return static_cast<std::uint64_t>(
        std::count_if(values.begin(), values.end(), predicate));
  };
  return {count(scores.electron,
                [&](double score) { return score > threshold; }),
          count(scores.pim,
                [&](double score) { return score <= threshold; }),
          count(scores.pi0,
                [&](double score) { return score <= threshold; }),
          static_cast<std::uint64_t>(scores.electron.size())};
}

double wilsonUpper(std::uint64_t selected, std::uint64_t denominator) {
  constexpr double z = 1.959963984540054;
  const double n = static_cast<double>(denominator);
  const double p = static_cast<double>(selected) / n;
  const double scale = 1.0 + z * z / n;
  const double centre = (p + z * z / (2.0 * n)) / scale;
  const double half = z * std::sqrt(
      p * (1.0 - p) / n + z * z / (4.0 * n * n)) / scale;
  return centre + half;
}

double totalBackgroundToBeta(const Counts &counts, bool upper) {
  const std::uint64_t pimSurvive = counts.denominator - counts.pimReject;
  const std::uint64_t pi0Survive = counts.denominator - counts.pi0Reject;
  if (upper) {
    return kPi0ToBetaFactor * wilsonUpper(pi0Survive, counts.denominator) +
           kPimToBetaFactor * wilsonUpper(pimSurvive, counts.denominator);
  }
  return kPi0ToBetaFactor * static_cast<double>(pi0Survive) / counts.denominator +
         kPimToBetaFactor * static_cast<double>(pimSurvive) / counts.denominator;
}

double thresholdForTotalRequirement(const SampleScores &scores, bool upper) {
  std::vector<double> orderedPim = scores.pim;
  std::vector<double> orderedPi0 = scores.pi0;
  std::sort(orderedPim.begin(), orderedPim.end());
  std::sort(orderedPi0.begin(), orderedPi0.end());
  std::vector<double> candidates;
  candidates.reserve(orderedPim.size() + orderedPi0.size());
  candidates.insert(candidates.end(), orderedPim.begin(), orderedPim.end());
  candidates.insert(candidates.end(), orderedPi0.begin(), orderedPi0.end());
  std::sort(candidates.begin(), candidates.end());
  candidates.erase(std::unique(candidates.begin(), candidates.end()), candidates.end());
  for (double threshold : candidates) {
    const auto pimAccepted = static_cast<std::uint64_t>(
        orderedPim.end() - std::upper_bound(
            orderedPim.begin(), orderedPim.end(), threshold));
    const auto pi0Accepted = static_cast<std::uint64_t>(
        orderedPi0.end() - std::upper_bound(
            orderedPi0.begin(), orderedPi0.end(), threshold));
    const auto n = static_cast<std::uint64_t>(orderedPim.size());
    const double total = upper
        ? kPi0ToBetaFactor * wilsonUpper(pi0Accepted, n) +
              kPimToBetaFactor * wilsonUpper(pimAccepted, n)
        : kPi0ToBetaFactor * static_cast<double>(pi0Accepted) / n +
              kPimToBetaFactor * static_cast<double>(pimAccepted) / n;
    if (total < kTotalRequirement) {
      return threshold;
    }
  }
  return candidates.back();
}

Counts baselineCounts(const fs::path &directory,
                      const GeometrySpec &geometry, int split) {
  const auto countSpecies = [&](const std::string &species, bool keep) {
    std::uint64_t selected = 0;
    std::uint64_t total = 0;
    const fs::path path = directory /
        (geometry.prefix + "_balanced_pc3_h60_s7232026_" + species + ".bgo2");
    forEachRow(path, geometry,
               [&](std::int64_t eventId, const std::vector<float> &row) {
      if ((eventId & 3) != split) return;
      ++total;
      const bool accepted =
          static_cast<int>(std::llround(row[4])) == 1 && row[kPcUpColumn] < 0.2f;
      if ((keep && accepted) || (!keep && !accepted)) ++selected;
    });
    if (total != 25000) throw std::runtime_error("unexpected split size");
    return selected;
  };
  return {countSpecies("e", true), countSpecies("pim", false),
          countSpecies("pi0", false), 25000};
}

void writeCounts(std::ostream &output, const Counts &counts) {
  output << "{\"denominator\":" << counts.denominator
         << ",\"electron_keep_count\":" << counts.electronKeep
         << ",\"electron_keep\":"
         << static_cast<double>(counts.electronKeep) / counts.denominator
         << ",\"pim_reject_count\":" << counts.pimReject
         << ",\"pim_reject\":"
         << static_cast<double>(counts.pimReject) / counts.denominator
         << ",\"pi0_reject_count\":" << counts.pi0Reject
         << ",\"pi0_reject\":"
         << static_cast<double>(counts.pi0Reject) / counts.denominator << '}';
}

void writeOperatingPoint(std::ostream &output, const SampleScores &scores,
                         double threshold) {
  const Counts counts = scoreCounts(scores, threshold);
  output << "{\"threshold\":" << threshold << ",\"metrics\":";
  writeCounts(output, counts);
  output << ",\"pi_total_to_beta\":"
         << totalBackgroundToBeta(counts, false)
         << ",\"pi_total_to_beta_mc_stat_upper95\":"
         << totalBackgroundToBeta(counts, true) << '}';
}

void writeScores(const fs::path &path, const SampleScores &scores,
                 double threshold) {
  std::ofstream output(path, std::ios::binary);
  if (!output) throw std::runtime_error("cannot create " + path.string());
  const std::array<char, 4> magic = {'B', 'D', 'T', '2'};
  const std::uint32_t version = 2;
  const std::array<std::uint32_t, 3> counts = {
      static_cast<std::uint32_t>(scores.electron.size()),
      static_cast<std::uint32_t>(scores.pim.size()),
      static_cast<std::uint32_t>(scores.pi0.size())};
  output.write(magic.data(), magic.size());
  output.write(reinterpret_cast<const char *>(&version), sizeof(version));
  output.write(reinterpret_cast<const char *>(counts.data()), sizeof(counts));
  output.write(reinterpret_cast<const char *>(&threshold), sizeof(threshold));
  for (const auto *values : {&scores.electron, &scores.pim, &scores.pi0}) {
    output.write(reinterpret_cast<const char *>(values->data()),
                 static_cast<std::streamsize>(values->size() * sizeof(double)));
  }
  if (!output) throw std::runtime_error("failed writing " + path.string());
}

}  // namespace

int main(int argc, char **argv) {
  try {
    if (argc != 4 && argc != 5) {
      std::cerr << "usage: balanced_pc3_pim_bdt_study PROJECT INPUT_BGO2_DIR "
                   "SCRATCH_WORK_DIR\n";
      return 2;
    }
    const fs::path project = fs::path(argv[1]);
    const fs::path inputs = fs::path(argv[2]);
    const fs::path work = fs::path(argv[3]);
    const std::string runMode = argc == 5 ? std::string(argv[4]) : "--train-all";
    if (runMode != "--train-all" && runMode != "--reuse" &&
        runMode != "--train-missing") {
      throw std::runtime_error("run mode must be --train-all, --reuse, or --train-missing");
    }
    const fs::path models = project / "analysis/models";
    const fs::path result = project /
        "analysis/results/bgoegg_v1/balanced_pc3_pim_bdt_s7232026.json";
    fs::create_directories(work);
    fs::create_directories(models);
    fs::create_directories(result.parent_path());
    TMVA::Tools::Instance();

    std::ofstream json(result);
    if (!json) throw std::runtime_error("cannot create " + result.string());
    json << std::setprecision(12);
    json << "{\n  \"status\":\"within-seed untouched-test study\",\n"
         << "  \"tag\":\"balanced_pc3_pim_bdt_s7232026\",\n"
         << "  \"split\":{\"training\":\"eventID mod4 0 or 1, 50000/species\","
         << "\"validation\":\"eventID mod4 2, 25000/species\","
         << "\"test\":\"eventID mod4 3, 25000/species\"},\n"
         << "  \"model\":{\"method\":\"" << kMethod
         << "\",\"options\":\"" << kOptions << "\"},\n"
         << "  \"cut_rule\":\"accept electron if BDT score > threshold; "
            "BR4 thresholds fixed on validation to maximize electron keep subject to "
            "pi total / beta < 0.04\",\n"
         << "  \"branching_requirement\":{\"beta_br\":0.00048,"
            "\"pi0_factor\":3.33333333333333,\"pim_factor\":0.0375,"
            "\"formula\":\"3.33333333333333*pi0_survival + 0.0375*pim_survival\","
            "\"maximum\":0.04,\"source\":\"Hong thesis Table 2.2 and Kamada thesis pp.37-38\"},\n"
         << "  \"truth_usage\":\"species label for training only; no truth variable input\",\n"
         << "  \"geometries\":{\n";

    for (std::size_t geometryIndex = 0; geometryIndex < kGeometries.size();
         ++geometryIndex) {
      const auto &geometry = kGeometries[geometryIndex];
      const Counts validationBaseline = baselineCounts(inputs, geometry, 2);
      const Counts testBaseline = baselineCounts(inputs, geometry, 3);
      if (geometryIndex) json << ",\n";
      json << "    \"" << geometry.name << "\":{\n      \"baseline_validation\":";
      writeCounts(json, validationBaseline);
      json << ",\n      \"baseline_test\":";
      writeCounts(json, testBaseline);
      json << ",\n      \"models\":{\n";

      for (int modelIndex = 0; modelIndex < 3; ++modelIndex) {
        const bool includePc = modelIndex != 0;
        const bool weightedTotal = modelIndex == 2;
        const std::string featureSet = modelIndex == 0
            ? "bgo_only" : (modelIndex == 1 ? "bgo_pc" : "bgo_pc_brweighted");
        const std::string factoryPrefix = weightedTotal ? "PITOTAL_" : "PIM_";
        const std::size_t nFeature = includePc ? 24 : 23;
        const fs::path modelWork = work / (geometry.prefix + "_" + featureSet);
        fs::create_directories(modelWork);
        const auto inputPath = [&](const std::string &species) {
          return inputs / (geometry.prefix +
              "_balanced_pc3_h60_s7232026_" + species + ".bgo2");
        };
        const fs::path trainingRoot = modelWork / "tmva_training.root";
        const fs::path generatedWeight = modelWork / "dataset/weights" /
            (factoryPrefix + geometry.prefix + "_" + featureSet + "_" + kMethod +
             ".weights.xml");
        const fs::path savedWeight = models /
            ("balanced_pc3_" + geometry.prefix + "_" + featureSet +
             "_bdtg_s7232026.weights.xml");
        const bool trainThis = runMode == "--train-all" ||
            (runMode == "--train-missing" && !fs::exists(savedWeight));
        if (runMode == "--reuse" && !fs::exists(savedWeight)) {
          throw std::runtime_error("missing saved weight " + savedWeight.string());
        }
        if (trainThis) {
          OwnedTree trainElectron("train_e", nFeature);
          OwnedTree trainPim("train_pim", nFeature);
          OwnedTree trainPi0("train_pi0", nFeature);
          fillTree(inputPath("e"), geometry, includePc, trainElectron);
          fillTree(inputPath("pim"), geometry, includePc, trainPim);
          if (weightedTotal) {
            fillTree(inputPath("pi0"), geometry, includePc, trainPi0);
          }
          if (trainElectron.tree->GetEntries() != 50000 ||
              trainPim.tree->GetEntries() != 50000 ||
              (weightedTotal && trainPi0.tree->GetEntries() != 50000)) {
            throw std::runtime_error("unexpected training split size");
          }
          gRandom->SetSeed(
              7232026 + static_cast<unsigned>(100 * geometryIndex + modelIndex));
          const fs::path originalDirectory = fs::current_path();
          fs::current_path(modelWork);
          {
            auto output = std::unique_ptr<TFile>(
                TFile::Open("tmva_training.root", "RECREATE"));
            if (!output || output->IsZombie()) {
              throw std::runtime_error("cannot create " + trainingRoot.string());
            }
            TMVA::Factory factory(
                (factoryPrefix + geometry.prefix + "_" + featureSet).c_str(),
                output.get(),
                "!V:!Silent:!Color:!DrawProgressBar:AnalysisType=Classification");
            TMVA::DataLoader loader("dataset");
            for (std::size_t index = 0; index < nFeature; ++index) {
              loader.AddVariable(kFeatures[index].c_str(), 'F');
            }
            loader.AddSignalTree(
                trainElectron.tree.get(), 1.0, TMVA::Types::kTraining);
            if (weightedTotal) {
              loader.AddBackgroundTree(
                  trainPim.tree.get(), kPimToBetaFactor,
                  TMVA::Types::kTraining);
              loader.AddBackgroundTree(
                  trainPi0.tree.get(), kPi0ToBetaFactor,
                  TMVA::Types::kTraining);
            } else {
              loader.AddBackgroundTree(
                  trainPim.tree.get(), 1.0, TMVA::Types::kTraining);
            }
            loader.PrepareTrainingAndTestTree(
                "", "", "NormMode=EqualNumEvents:!V");
            factory.BookMethod(&loader, TMVA::Types::kBDT, kMethod, kOptions);
            factory.TrainAllMethods();
            output->Close();
          }
          fs::current_path(originalDirectory);
          fs::copy_file(generatedWeight, savedWeight,
                        fs::copy_options::overwrite_existing);
        }

        TMVA::Reader reader("!Color:!Silent");
        std::array<float, kMaximumFeatureCount> values{};
        for (std::size_t index = 0; index < nFeature; ++index) {
          reader.AddVariable(kFeatures[index].c_str(), &values[index]);
        }
        reader.BookMVA(kMethod, savedWeight.c_str());
        const SampleScores train = evaluateSample(
            inputs, geometry, includePc, 0, reader, values);
        const SampleScores validation = evaluateSample(
            inputs, geometry, includePc, 2, reader, values);
        const SampleScores test = evaluateSample(
            inputs, geometry, includePc, 3, reader, values);
        const double targetKeep = static_cast<double>(validationBaseline.electronKeep) /
                                  validationBaseline.denominator;
        const double threshold =
            thresholdForElectronKeep(validation.electron, targetKeep);
        const Counts validationCounts = scoreCounts(validation, threshold);
        const Counts testCounts = scoreCounts(test, threshold);
        const double br4CentralThreshold =
            thresholdForTotalRequirement(validation, false);
        const double br4GuardedThreshold =
            thresholdForTotalRequirement(validation, true);
        const fs::path scorePath =
            work / (geometry.prefix + "_" + featureSet + "_test.bdt2");
        writeScores(scorePath, test, threshold);

        if (modelIndex) json << ",\n";
        json << "        \"" << featureSet << "\":{\"n_features\":"
             << nFeature << ",\"training_background\":\""
             << (weightedTotal
                     ? "pi0 and pim weighted by post-TH/TLC/delta-z BR factors"
                     : "pim only")
             << "\",\"features\":[";
        for (std::size_t index = 0; index < nFeature; ++index) {
          if (index) json << ',';
          json << '\"' << kFeatures[index] << '\"';
        }
        json << "],\"threshold\":" << threshold
             << ",\"weight_file\":\"" << savedWeight.string()
             << "\",\"score_file\":\"" << scorePath.string()
             << "\",\"auc_train_mod4_0\":"
             << auc(train.electron, train.pim)
             << ",\"auc_validation\":"
             << auc(validation.electron, validation.pim)
             << ",\"auc_test\":" << auc(test.electron, test.pim)
             << ",\"validation\":";
        writeCounts(json, validationCounts);
        json << ",\"test\":";
        writeCounts(json, testCounts);
        json << ",\"br4_central\":{\"validation\":";
        writeOperatingPoint(json, validation, br4CentralThreshold);
        json << ",\"test\":";
        writeOperatingPoint(json, test, br4CentralThreshold);
        json << "},\"br4_guarded\":{\"validation\":";
        writeOperatingPoint(json, validation, br4GuardedThreshold);
        json << ",\"test\":";
        writeOperatingPoint(json, test, br4GuardedThreshold);
        json << '}';
        json << '}';

        std::cout << geometry.name << ' ' << featureSet
                  << " test: e keep="
                  << 100.0 * testCounts.electronKeep / testCounts.denominator
                  << "% pim reject="
                  << 100.0 * testCounts.pimReject / testCounts.denominator
                  << "% pi0 reject="
                  << 100.0 * testCounts.pi0Reject / testCounts.denominator
                  << "% AUC=" << auc(test.electron, test.pim) << '\n';
        const Counts br4Test = scoreCounts(test, br4CentralThreshold);
        std::cout << geometry.name << ' ' << featureSet
                  << " BR4 test: e keep="
                  << 100.0 * br4Test.electronKeep / br4Test.denominator
                  << "% total/beta="
                  << 100.0 * totalBackgroundToBeta(br4Test, false)
                  << "% upper95="
                  << 100.0 * totalBackgroundToBeta(br4Test, true) << "%\n";
      }
      json << "\n      }\n    }";
    }
    json << "\n  }\n}\n";
    std::cout << result << '\n';
    return 0;
  } catch (const std::exception &error) {
    std::cerr << "ERROR: " << error.what() << '\n';
    return 1;
  }
}
