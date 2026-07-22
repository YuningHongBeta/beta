#include <TMVA/Reader.h>

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
constexpr std::size_t kFeatureCount = 24;
constexpr std::size_t kPcUpColumn = 64;
constexpr double kThreshold = 0.953974177251;
constexpr double kPi0ToBeta = 0.2 * 0.1 * 0.1 * 0.8 / 0.00048;
constexpr double kPimToBeta = 0.4 * 0.001 * 0.045 / 0.00048;
constexpr const char *kMethod = "BDTG_D3_M1";

struct Bgo2Header {
  std::uint32_t magic, version, nrow, ncol;
  std::int32_t nLayer, nSector, segmentationMode, physicsFlag, nSegTH, nSegTLC;
  double thetaMinDeg, thetaMaxDeg, thresholdMeV;
};
static_assert(sizeof(Bgo2Header) == 64, "unexpected BGO2 header layout");

const std::array<std::string, kFeatureCount> kFeatures = {
    "maxE", "sumE", "nHit", "nCl4", "nCl8", "cl1Sum4", "cl2Sum4",
    "cl1Size4", "cl1MaxFrac4", "cl1RmsDeg4", "allRmsDeg", "cl2OverCl1",
    "isolatedEFrac4", "seed15NCl4", "seed15Cl1Sum4", "seed15Cl1Size4",
    "e2OverE1", "nHit2T", "sum2T", "local3x3Frac", "meanHitE",
    "thetaSpan", "phiSpan", "pcUpE_MeV"};

struct Selected {
  std::uint64_t electron{}, pim{}, pi0{}, denominator{100000};
};

template <class Callback>
void forEachRow(const fs::path &path, int layers, int sectors,
                int segmentation, Callback callback) {
  std::ifstream input(path, std::ios::binary);
  if (!input) throw std::runtime_error("cannot open " + path.string());
  Bgo2Header header{};
  input.read(reinterpret_cast<char *>(&header), sizeof(header));
  if (!input || header.magic != kMagic || header.version != kVersion ||
      header.nrow != 100000 || header.ncol <= kPcUpColumn ||
      header.nLayer != layers || header.nSector != sectors ||
      header.segmentationMode != segmentation || header.physicsFlag != 4 ||
      std::abs(header.thresholdMeV - 3.0) > 1e-12) {
    throw std::runtime_error("incompatible BGO2 header in " + path.string());
  }
  std::vector<float> row(header.ncol);
  for (std::uint32_t entry = 0; entry < header.nrow; ++entry) {
    input.read(reinterpret_cast<char *>(row.data()),
               static_cast<std::streamsize>(row.size() * sizeof(float)));
    if (!input) throw std::runtime_error("truncated row in " + path.string());
    callback(row);
  }
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

double totalToBeta(const Selected &selected, bool upper) {
  if (upper) {
    return kPi0ToBeta * wilsonUpper(selected.pi0, selected.denominator) +
           kPimToBeta * wilsonUpper(selected.pim, selected.denominator);
  }
  return kPi0ToBeta * static_cast<double>(selected.pi0) / selected.denominator +
         kPimToBeta * static_cast<double>(selected.pim) / selected.denominator;
}

Selected evaluateBgoc(const fs::path &directory) {
  Selected result;
  for (const std::string species : {"e", "pim", "pi0"}) {
    std::uint64_t count = 0;
    forEachRow(directory / ("bgoc_" + species + ".bgo2"), 15, 15, 0,
               [&](const std::vector<float> &row) {
      if (static_cast<int>(std::llround(row[4])) == 1 &&
          row[kPcUpColumn] < 0.2f) ++count;
    });
    if (species == "e") result.electron = count;
    if (species == "pim") result.pim = count;
    if (species == "pi0") result.pi0 = count;
  }
  return result;
}

Selected evaluateBgoegg(const fs::path &directory, const fs::path &weightPath) {
  TMVA::Reader reader("!Color:!Silent");
  std::array<float, kFeatureCount> values{};
  for (std::size_t index = 0; index < kFeatureCount; ++index)
    reader.AddVariable(kFeatures[index].c_str(), &values[index]);
  reader.BookMVA(kMethod, weightPath.c_str());
  Selected result;
  for (const std::string species : {"e", "pim", "pi0"}) {
    std::uint64_t count = 0;
    forEachRow(directory / ("bgoegg31_" + species + ".bgo2"), 31, 60, 2,
               [&](const std::vector<float> &row) {
      for (std::size_t index = 0; index < kFeatureCount - 1; ++index)
        values[index] = row[index + 1];
      values[kFeatureCount - 1] = row[kPcUpColumn];
      if (reader.EvaluateMVA(kMethod) > kThreshold) ++count;
    });
    if (species == "e") result.electron = count;
    if (species == "pim") result.pim = count;
    if (species == "pi0") result.pi0 = count;
  }
  return result;
}

void writeResult(std::ostream &out, const Selected &s) {
  out << "{\"denominator_per_species\":" << s.denominator
      << ",\"electron_keep_count\":" << s.electron
      << ",\"electron_keep\":" << static_cast<double>(s.electron) / s.denominator
      << ",\"pim_reject_count\":" << s.denominator - s.pim
      << ",\"pim_reject\":" << 1.0 - static_cast<double>(s.pim) / s.denominator
      << ",\"pi0_reject_count\":" << s.denominator - s.pi0
      << ",\"pi0_reject\":" << 1.0 - static_cast<double>(s.pi0) / s.denominator
      << ",\"pi_total_to_beta\":" << totalToBeta(s, false)
      << ",\"pi_total_to_beta_mc_stat_upper95\":" << totalToBeta(s, true) << '}';
}

void printResult(const std::string &geometry, const Selected &s) {
  std::cout << geometry << ": e=" << 100.0 * s.electron / s.denominator
            << "% pim reject=" << 100.0 * (1.0 - static_cast<double>(s.pim) / s.denominator)
            << "% pi0 reject=" << 100.0 * (1.0 - static_cast<double>(s.pi0) / s.denominator)
            << "% total/beta=" << 100.0 * totalToBeta(s, false)
            << "% upper95=" << 100.0 * totalToBeta(s, true) << "%\n";
}
}  // namespace

int main(int argc, char **argv) {
  try {
    if (argc != 4) {
      std::cerr << "usage: balanced_pc3_ml_confirm PROJECT INPUT_BGO2_DIR OUTPUT_JSON\n";
      return 2;
    }
    const fs::path project = argv[1];
    const fs::path inputs = argv[2];
    const fs::path outputPath = argv[3];
    const fs::path weightPath = project /
        "analysis/models/balanced_pc3_bgoegg31_bgo_pc_brweighted_bdtg_s7232026.weights.xml";
    const Selected bgoc = evaluateBgoc(inputs);
    const Selected bgoegg = evaluateBgoegg(inputs, weightPath);
    fs::create_directories(outputPath.parent_path());
    std::ofstream output(outputPath);
    if (!output) throw std::runtime_error("cannot create " + outputPath.string());
    output << std::setprecision(12)
           << "{\n  \"status\":\"independent-seed frozen-cut confirmation\",\n"
           << "  \"tag\":\"balanced_pc3_ml_confirm_s7242026\",\n"
           << "  \"development_seed\":7232026,\n"
           << "  \"confirmation_seed\":7242026,\n"
           << "  \"events_per_species\":100000,\n"
           << "  \"branching_requirement\":{\"beta_br\":0.00048,"
              "\"pi0_factor\":3.33333333333333,\"pim_factor\":0.0375,"
              "\"formula\":\"3.33333333333333*pi0_survival + 0.0375*pim_survival\","
              "\"maximum\":0.04},\n"
           << "  \"selection\":{\"BGOC\":\"T3 nCl4 == 1 and pcUpE_MeV < 0.2\","
              "\"BGOegg31\":\"BR-weighted BDTG; score > frozen guarded-validation threshold\"},\n"
           << "  \"bgoegg31_frozen_threshold\":" << kThreshold << ",\n"
           << "  \"results\":{\"BGOC\":";
    writeResult(output, bgoc);
    output << ",\"BGOegg31\":";
    writeResult(output, bgoegg);
    output << "}\n}\n";
    printResult("BGOC", bgoc);
    printResult("BGOegg31", bgoegg);
    return 0;
  } catch (const std::exception &error) {
    std::cerr << "ERROR: " << error.what() << '\n';
    return 1;
  }
}
