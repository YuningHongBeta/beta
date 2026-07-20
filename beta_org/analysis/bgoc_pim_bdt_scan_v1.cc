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
#include <utility>
#include <vector>

namespace fs = std::filesystem;

namespace {

constexpr std::uint32_t kMagic = 0x42474f32U;
constexpr std::uint32_t kVersion = 2U;
constexpr double kNominalTarget = 0.9000;
constexpr double kGuardedTarget = 0.9050;
constexpr double kFrozenTarget = 0.9040;
constexpr double kFrozenThreshold = 0.480194806111;
constexpr const char *kFrozenMethod = "BDTG_D3_M1";
const std::array<double, 8> kTargetScan = {
    0.900, 0.901, 0.902, 0.903, 0.904, 0.905, 0.906, 0.908};

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

struct MethodSpec {
  std::string name;
  std::string options;
};

const std::vector<MethodSpec> kMethods = {
    {"BDTG_D2",
     "!H:!V:NTrees=500:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.05:"
     "UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=30:MaxDepth=2"},
    {"BDTG_D3",
     "!H:!V:NTrees=500:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.05:"
     "UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=30:MaxDepth=3"},
    {"BDTG_D4",
     "!H:!V:NTrees=500:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.05:"
     "UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=30:MaxDepth=4"},
    {"BDTG_D3_M1",
     "!H:!V:NTrees=800:MinNodeSize=1.0%:BoostType=Grad:Shrinkage=0.03:"
     "UseBaggedBoost:BaggedSampleFraction=0.6:nCuts=40:MaxDepth=3"},
    {"BDT_Ada_D2",
     "!H:!V:NTrees=500:MinNodeSize=2.5%:BoostType=AdaBoost:AdaBoostBeta=0.5:"
     "UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=30:MaxDepth=2"},
    {"BDT_Ada_D3",
     "!H:!V:NTrees=500:MinNodeSize=2.5%:BoostType=AdaBoost:AdaBoostBeta=0.5:"
     "UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=30:MaxDepth=3"},
};

struct OwnedTree {
  std::unique_ptr<TTree> tree;
  std::array<float, 23> values{};

  explicit OwnedTree(const std::string &name)
      : tree(std::make_unique<TTree>(name.c_str(), name.c_str())) {
    for (std::size_t i = 0; i < kFeatures.size(); ++i) {
      tree->Branch(kFeatures[i].c_str(), &values[i],
                   (kFeatures[i] + "/F").c_str());
    }
  }
};

struct Scores {
  std::vector<double> electron;
  std::vector<double> pim;
  std::vector<double> pi0;
};

struct Metrics {
  double electronKeep{};
  double pimReject{};
  double pi0Reject{};
  double aucElectronPim{};
};

template <class Callback>
void forEachBgo2Row(const fs::path &path, Callback callback) {
  std::ifstream input(path, std::ios::binary);
  if (!input) {
    throw std::runtime_error("cannot open " + path.string());
  }
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
    if (!input) {
      throw std::runtime_error("truncated BGO2 row in " + path.string());
    }
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

void fillTree(const fs::path &path, int splitValue, OwnedTree &output) {
  forEachBgo2Row(path, [&](std::int64_t eventId, const std::vector<float> &row) {
    if ((eventId & 3) != splitValue) {
      return;
    }
    for (std::size_t i = 0; i < kFeatures.size(); ++i) {
      output.values[i] = row[i + 1];
    }
    output.tree->Fill();
  });
}

double thresholdForReject(const std::vector<double> &background, double target) {
  if (background.empty()) {
    throw std::runtime_error("empty background score vector");
  }
  std::vector<double> ordered = background;
  std::sort(ordered.begin(), ordered.end());
  const auto index = std::min<std::size_t>(
      ordered.size() - 1,
      static_cast<std::size_t>(std::ceil(target * ordered.size())) - 1);
  return ordered[index];
}

double auc(const std::vector<double> &electron,
           const std::vector<double> &background) {
  struct Item {
    double score;
    bool isBackground;
  };
  std::vector<Item> pooled;
  pooled.reserve(electron.size() + background.size());
  for (double value : electron) pooled.push_back({value, false});
  for (double value : background) pooled.push_back({value, true});
  std::sort(pooled.begin(), pooled.end(),
            [](const Item &a, const Item &b) { return a.score < b.score; });
  double backgroundRankSum = 0.0;
  std::size_t begin = 0;
  while (begin < pooled.size()) {
    std::size_t end = begin + 1;
    while (end < pooled.size() && pooled[end].score == pooled[begin].score) ++end;
    const double averageRank = 0.5 * (static_cast<double>(begin + 1) +
                                      static_cast<double>(end));
    for (std::size_t i = begin; i < end; ++i) {
      if (pooled[i].isBackground) backgroundRankSum += averageRank;
    }
    begin = end;
  }
  const double nb = static_cast<double>(background.size());
  const double ne = static_cast<double>(electron.size());
  // TMVA score is electron-like when large. Return P(e score > pim score).
  const double uBackground = backgroundRankSum - nb * (nb + 1.0) / 2.0;
  return 1.0 - uBackground / (ne * nb);
}

Metrics metrics(const Scores &scores, double threshold) {
  const auto fraction = [](const std::vector<double> &values, auto predicate) {
    return static_cast<double>(std::count_if(values.begin(), values.end(), predicate)) /
           static_cast<double>(values.size());
  };
  return {
      fraction(scores.electron, [&](double value) { return value > threshold; }),
      fraction(scores.pim, [&](double value) { return value <= threshold; }),
      fraction(scores.pi0, [&](double value) { return value <= threshold; }),
      auc(scores.electron, scores.pim),
  };
}

std::vector<double> evaluateFile(const fs::path &path, TMVA::Reader &reader,
                                 const std::string &method,
                                 std::array<float, 23> &values,
                                 int splitValue) {
  std::vector<double> result;
  forEachBgo2Row(path, [&](std::int64_t eventId, const std::vector<float> &row) {
    if (splitValue >= 0 && (eventId & 3) != splitValue) return;
    for (std::size_t i = 0; i < kFeatures.size(); ++i) values[i] = row[i + 1];
    result.push_back(reader.EvaluateMVA(method.c_str()));
  });
  return result;
}

Scores evaluateSample(const fs::path &directory, TMVA::Reader &reader,
                      const std::string &method, std::array<float, 23> &values,
                      int splitValue) {
  return {
      evaluateFile(directory / "e.bgo2", reader, method, values, splitValue),
      evaluateFile(directory / "pim.bgo2", reader, method, values, splitValue),
      evaluateFile(directory / "pi0.bgo2", reader, method, values, splitValue),
  };
}

void writeMetric(std::ostream &out, const Metrics &metric) {
  out << "{\"electron_keep\":" << metric.electronKeep
      << ",\"pim_reject\":" << metric.pimReject
      << ",\"pi0_reject\":" << metric.pi0Reject
      << ",\"auc_e_pim\":" << metric.aucElectronPim << "}";
}

}  // namespace

int main(int argc, char **argv) {
  try {
    const fs::path project = argc > 1 ? fs::path(argv[1]) : fs::current_path();
    const bool reuseWeights = argc > 2 && std::string(argv[2]) == "--reuse";
    const fs::path training =
        project / "tmp/bgoc_u15x15_feasibility_v1/bgo2/none/training_T3";
    const fs::path dev730 =
        project / "tmp/bgoc_u15x15_feasibility_v1/bgo2/none/confirmation_T3";
    const fs::path dev930 = project / "tmp/bgoc_qda5_confirm_s9302026";
    const fs::path work = project / "tmp/bgoc_pim_bdt_scan_v1";
    const fs::path resultPath =
        project / "analysis/results/bgoegg_v1/bgoc_pim_bdt_scan_v1.json";
    fs::create_directories(work);
    fs::create_directories(resultPath.parent_path());

    OwnedTree trainElectron("train_e");
    OwnedTree trainPim("train_pim");
    OwnedTree validationElectron("validation_e");
    OwnedTree validationPim("validation_pim");
    fillTree(training / "e.bgo2", 0, trainElectron);
    fillTree(training / "pim.bgo2", 0, trainPim);
    fillTree(training / "e.bgo2", 2, validationElectron);
    fillTree(training / "pim.bgo2", 2, validationPim);
    if (trainElectron.tree->GetEntries() != 25000 ||
        trainPim.tree->GetEntries() != 25000 ||
        validationElectron.tree->GetEntries() != 25000 ||
        validationPim.tree->GetEntries() != 25000) {
      throw std::runtime_error("unexpected seed630 split counts");
    }

    TMVA::Tools::Instance();
    const fs::path trainingRoot = work / "tmva_training.root";
    if (!reuseWeights) {
      gRandom->SetSeed(6302026);
      auto output = std::unique_ptr<TFile>(TFile::Open(trainingRoot.c_str(), "RECREATE"));
      if (!output || output->IsZombie()) {
        throw std::runtime_error("cannot create TMVA output ROOT");
      }
      TMVA::Factory factory(
          "BGOC_PIM", output.get(),
          "!V:!Silent:Color:DrawProgressBar:AnalysisType=Classification");
      TMVA::DataLoader loader((work / "dataset").c_str());
      for (const auto &feature : kFeatures) loader.AddVariable(feature, 'F');
      loader.AddSignalTree(trainElectron.tree.get(), 1.0, TMVA::Types::kTraining);
      loader.AddBackgroundTree(trainPim.tree.get(), 1.0, TMVA::Types::kTraining);
      loader.AddSignalTree(validationElectron.tree.get(), 1.0,
                           TMVA::Types::kTesting);
      loader.AddBackgroundTree(validationPim.tree.get(), 1.0,
                               TMVA::Types::kTesting);
      loader.PrepareTrainingAndTestTree(
          "", "", "NormMode=EqualNumEvents:!V");
      for (const auto &method : kMethods) {
        factory.BookMethod(&loader, TMVA::Types::kBDT, method.name.c_str(),
                           method.options.c_str());
      }
      factory.TrainAllMethods();
      factory.TestAllMethods();
      factory.EvaluateAllMethods();
      output->Close();
    }

    std::ofstream json(resultPath);
    if (!json) throw std::runtime_error("cannot create result JSON");
    json << std::setprecision(12);
    json << "{\n  \"status\":\"development_scan\",\n"
         << "  \"analysis\":\"bgoc_pim_bdt_scan_v1\",\n"
         << "  \"objective\":\"pim reject > 90% and electron keep > 60%\",\n"
         << "  \"training\":\"seed6302026 eventID mod4=0, 25000/species\",\n"
         << "  \"cut_selection\":\"seed6302026 eventID mod4=2, 25000/species\",\n"
         << "  \"development_audits\":[\"seed7302026\",\"seed9302026\"],\n"
         << "  \"pi0_role\":\"guard only; not used in BDT training\",\n"
         << "  \"frozen_selection\":{\"method\":\"" << kFrozenMethod
         << "\",\"selection_target\":" << kFrozenTarget
         << ",\"threshold\":" << kFrozenThreshold
         << ",\"cut_rule\":\"accept electron if score > threshold\","
         << "\"basis\":\"passes strict electron keep > 60% and pim reject > 90% on cut-selection, seed730 audit, and seed930 audit\"},\n"
         << "  \"features\":[";
    for (std::size_t i = 0; i < kFeatures.size(); ++i) {
      if (i) json << ',';
      json << '\"' << kFeatures[i] << '\"';
    }
    json << "],\n  \"methods\":[\n";

    for (std::size_t index = 0; index < kMethods.size(); ++index) {
      const auto &method = kMethods[index];
      TMVA::Reader reader("!Color:!Silent");
      std::array<float, 23> values{};
      for (std::size_t i = 0; i < kFeatures.size(); ++i) {
        reader.AddVariable(kFeatures[i].c_str(), &values[i]);
      }
      const fs::path weight = work / "dataset/weights" /
                              ("BGOC_PIM_" + method.name + ".weights.xml");
      reader.BookMVA(method.name.c_str(), weight.c_str());
      const Scores validation =
          evaluateSample(training, reader, method.name, values, 2);
      const Scores audit730 =
          evaluateSample(dev730, reader, method.name, values, -1);
      const Scores audit930 =
          evaluateSample(dev930, reader, method.name, values, -1);
      const double nominalThreshold =
          thresholdForReject(validation.pim, kNominalTarget);
      const double guardedThreshold =
          thresholdForReject(validation.pim, kGuardedTarget);

      if (index) json << ",\n";
      json << "    {\"name\":\"" << method.name << "\",\"options\":\""
           << method.options << "\",\"nominal_threshold\":"
           << nominalThreshold << ",\"guarded_threshold\":" << guardedThreshold
           << ",\"nominal\":{\"validation\":";
      writeMetric(json, metrics(validation, nominalThreshold));
      json << ",\"audit730\":";
      writeMetric(json, metrics(audit730, nominalThreshold));
      json << ",\"audit930\":";
      writeMetric(json, metrics(audit930, nominalThreshold));
      json << "},\"guarded\":{\"validation\":";
      writeMetric(json, metrics(validation, guardedThreshold));
      json << ",\"audit730\":";
      writeMetric(json, metrics(audit730, guardedThreshold));
      json << ",\"audit930\":";
      writeMetric(json, metrics(audit930, guardedThreshold));
      json << "},\"operating_points\":[";
      for (std::size_t targetIndex = 0; targetIndex < kTargetScan.size();
           ++targetIndex) {
        const double target = kTargetScan[targetIndex];
        const double threshold = thresholdForReject(validation.pim, target);
        if (targetIndex) json << ',';
        json << "{\"selection_target\":" << target
             << ",\"threshold\":" << threshold << ",\"validation\":";
        writeMetric(json, metrics(validation, threshold));
        json << ",\"audit730\":";
        writeMetric(json, metrics(audit730, threshold));
        json << ",\"audit930\":";
        writeMetric(json, metrics(audit930, threshold));
        json << '}';
      }
      json << "]}";

      const Metrics check = metrics(audit930, guardedThreshold);
      std::cout << method.name << " guarded audit930: e keep="
                << 100.0 * check.electronKeep
                << "% pim reject=" << 100.0 * check.pimReject
                << "% pi0 reject=" << 100.0 * check.pi0Reject
                << "% AUC=" << check.aucElectronPim << '\n';
    }
    json << "\n  ]\n}\n";
    std::cout << resultPath << '\n';
    return 0;
  } catch (const std::exception &error) {
    std::cerr << "ERROR: " << error.what() << '\n';
    return 1;
  }
}
