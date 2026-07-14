#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr std::uint32_t kMagic = 0x42474f32U; // BGO2
constexpr std::uint32_t kVersion = 2U;
constexpr int kPatchThetaRadius = 3;
constexpr int kPatchPhiRadius = 5;

struct Meta {
  int nLayer = 0;
  int nSector = 0;
  int segmentationMode = -1;
  int physicsFlag = -1;
  int nSegTH = 0;
  int nSegTLC = 0;
  double thetaMinDeg = 0.0;
  double thetaMaxDeg = 0.0;
};

void requireBranch(TTree *tree, const char *name) {
  if (!tree || !tree->GetBranch(name))
    throw std::runtime_error(std::string("missing branch ") +
                             (tree ? tree->GetName() : "<null>") + "." + name);
}

Meta readMeta(TTree *tree) {
  if (!tree || tree->GetEntries() != 1)
    throw std::runtime_error("runmeta must contain exactly one entry");
  const char *required[] = {
      "nLayer", "nSector", "segmentationMode", "physicsFlag",
      "nSegTH", "nSegTLC", "thetaMin_deg", "thetaMax_deg",
  };
  for (const char *name : required) requireBranch(tree, name);
  Meta meta;
  tree->SetBranchStatus("*", 0);
  for (const char *name : required) tree->SetBranchStatus(name, 1);
  tree->SetBranchAddress("nLayer", &meta.nLayer);
  tree->SetBranchAddress("nSector", &meta.nSector);
  tree->SetBranchAddress("segmentationMode", &meta.segmentationMode);
  tree->SetBranchAddress("physicsFlag", &meta.physicsFlag);
  tree->SetBranchAddress("nSegTH", &meta.nSegTH);
  tree->SetBranchAddress("nSegTLC", &meta.nSegTLC);
  tree->SetBranchAddress("thetaMin_deg", &meta.thetaMinDeg);
  tree->SetBranchAddress("thetaMax_deg", &meta.thetaMaxDeg);
  tree->GetEntry(0);
  if (meta.nLayer <= 0 || meta.nSector <= 0 || meta.physicsFlag != 4)
    throw std::runtime_error("unsupported run metadata");
  return meta;
}

template <typename T> void writeValue(std::ofstream &out, T value) {
  out.write(reinterpret_cast<const char *>(&value), sizeof(value));
}

std::vector<std::string> featureNames(const Meta &meta) {
  std::vector<std::string> names = {"eventID"};
  for (int layer = 0; layer < meta.nLayer; ++layer)
    names.push_back("ringSumE_" + std::to_string(layer));
  for (int layer = 0; layer < meta.nLayer; ++layer)
    names.push_back("ringNHit_" + std::to_string(layer));
  for (int layer = 0; layer < meta.nLayer; ++layer)
    names.push_back("ringMaxE_" + std::to_string(layer));
  for (int relativePhi = 0; relativePhi < meta.nSector; ++relativePhi)
    names.push_back("phiRelSumE_" + std::to_string(relativePhi));
  for (int relativePhi = 0; relativePhi < meta.nSector; ++relativePhi)
    names.push_back("phiRelNHit_" + std::to_string(relativePhi));
  for (int deltaTheta = -kPatchThetaRadius;
       deltaTheta <= kPatchThetaRadius; ++deltaTheta) {
    for (int deltaPhi = -kPatchPhiRadius;
         deltaPhi <= kPatchPhiRadius; ++deltaPhi) {
      names.push_back("patchE_dt" + std::to_string(deltaTheta) + "_dp" +
                      std::to_string(deltaPhi));
    }
  }
  return names;
}

std::vector<float> makeFeatures(const Meta &meta, int eventID,
                                const std::vector<double> &raw,
                                double threshold) {
  const int nCell = meta.nLayer * meta.nSector;
  if (static_cast<int>(raw.size()) != nCell)
    throw std::runtime_error("calarr vector size mismatch");
  std::vector<double> energy(nCell, 0.0);
  int leadingCell = -1;
  double leadingEnergy = 0.0;
  for (int cell = 0; cell < nCell; ++cell) {
    if (raw[cell] < threshold) continue;
    energy[cell] = raw[cell];
    if (energy[cell] > leadingEnergy) {
      leadingEnergy = energy[cell];
      leadingCell = cell;
    }
  }
  const int leadingTheta = leadingCell >= 0 ? leadingCell / meta.nSector : 0;
  const int leadingPhi = leadingCell >= 0 ? leadingCell % meta.nSector : 0;
  const auto names = featureNames(meta);
  std::vector<float> output(names.size(), 0.0f);
  output[0] = static_cast<float>(eventID);
  std::size_t offset = 1;
  for (int layer = 0; layer < meta.nLayer; ++layer) {
    double sum = 0.0;
    for (int phi = 0; phi < meta.nSector; ++phi)
      sum += energy[layer * meta.nSector + phi];
    output[offset++] = static_cast<float>(sum);
  }
  for (int layer = 0; layer < meta.nLayer; ++layer) {
    int nHit = 0;
    for (int phi = 0; phi < meta.nSector; ++phi)
      nHit += energy[layer * meta.nSector + phi] > 0.0;
    output[offset++] = static_cast<float>(nHit);
  }
  for (int layer = 0; layer < meta.nLayer; ++layer) {
    double maximum = 0.0;
    for (int phi = 0; phi < meta.nSector; ++phi)
      maximum = std::max(maximum, energy[layer * meta.nSector + phi]);
    output[offset++] = static_cast<float>(maximum);
  }
  for (int relativePhi = 0; relativePhi < meta.nSector; ++relativePhi) {
    const int phi = (leadingPhi + relativePhi) % meta.nSector;
    double sum = 0.0;
    for (int layer = 0; layer < meta.nLayer; ++layer)
      sum += energy[layer * meta.nSector + phi];
    output[offset++] = static_cast<float>(sum);
  }
  for (int relativePhi = 0; relativePhi < meta.nSector; ++relativePhi) {
    const int phi = (leadingPhi + relativePhi) % meta.nSector;
    int nHit = 0;
    for (int layer = 0; layer < meta.nLayer; ++layer)
      nHit += energy[layer * meta.nSector + phi] > 0.0;
    output[offset++] = static_cast<float>(nHit);
  }
  for (int deltaTheta = -kPatchThetaRadius;
       deltaTheta <= kPatchThetaRadius; ++deltaTheta) {
    const int layer = leadingTheta + deltaTheta;
    for (int deltaPhi = -kPatchPhiRadius;
         deltaPhi <= kPatchPhiRadius; ++deltaPhi) {
      if (leadingCell < 0 || layer < 0 || layer >= meta.nLayer) {
        output[offset++] = 0.0f;
        continue;
      }
      const int phi =
          (leadingPhi + deltaPhi + meta.nSector) % meta.nSector;
      output[offset++] = static_cast<float>(energy[layer * meta.nSector + phi]);
    }
  }
  if (offset != output.size())
    throw std::runtime_error("internal feature-count mismatch");
  return output;
}

void writeManifest(const std::string &outputPath, const std::string &inputPath,
                   const Meta &meta, double threshold, std::uint32_t nrow,
                   const std::vector<std::string> &names) {
  std::ofstream out(outputPath + ".json");
  out << std::setprecision(15)
      << "{\n"
      << "  \"format\": \"BGO2\",\n"
      << "  \"version\": 2,\n"
      << "  \"analysis\": \"bgo_pattern_v1\",\n"
      << "  \"input\": \"" << inputPath << "\",\n"
      << "  \"nrow\": " << nrow << ",\n"
      << "  \"ncol\": " << names.size() << ",\n"
      << "  \"nLayer\": " << meta.nLayer << ",\n"
      << "  \"nSector\": " << meta.nSector << ",\n"
      << "  \"segmentationMode\": " << meta.segmentationMode << ",\n"
      << "  \"physicsFlag\": " << meta.physicsFlag << ",\n"
      << "  \"nSegTH\": " << meta.nSegTH << ",\n"
      << "  \"nSegTLC\": " << meta.nSegTLC << ",\n"
      << "  \"thetaMin_deg\": " << meta.thetaMinDeg << ",\n"
      << "  \"thetaMax_deg\": " << meta.thetaMaxDeg << ",\n"
      << "  \"threshold_MeV\": " << threshold << ",\n"
      << "  \"alignment\": \"leading-cell phi mapped to relative phi zero\",\n"
      << "  \"truthUsage\": \"none\",\n"
      << "  \"features\": [\n";
  for (std::size_t index = 0; index < names.size(); ++index)
    out << "    \"" << names[index] << "\""
        << (index + 1 == names.size() ? "\n" : ",\n");
  out << "  ]\n}\n";
}

void extract(const std::string &inputPath, const std::string &outputPath,
             double threshold) {
  if (!(threshold > 0.0)) throw std::runtime_error("threshold must be positive");
  TFile input(inputPath.c_str(), "READ");
  if (input.IsZombie()) throw std::runtime_error("cannot open " + inputPath);
  auto *calarr = dynamic_cast<TTree *>(input.Get("calarr"));
  auto *runmeta = dynamic_cast<TTree *>(input.Get("runmeta"));
  if (!calarr || !runmeta) throw std::runtime_error("calarr/runmeta missing");
  const Meta meta = readMeta(runmeta);
  requireBranch(calarr, "eventID");
  requireBranch(calarr, "dE_MeV");
  int eventID = -1;
  std::vector<double> *energy = nullptr;
  calarr->SetBranchStatus("*", 0);
  calarr->SetBranchStatus("eventID", 1);
  calarr->SetBranchStatus("dE_MeV", 1);
  calarr->SetBranchAddress("eventID", &eventID);
  calarr->SetBranchAddress("dE_MeV", &energy);
  const auto names = featureNames(meta);
  const auto nrow64 = calarr->GetEntries();
  if (nrow64 < 0 || nrow64 > 0xffffffffLL)
    throw std::runtime_error("unsupported row count");
  const auto nrow = static_cast<std::uint32_t>(nrow64);
  const auto ncol = static_cast<std::uint32_t>(names.size());
  std::ofstream output(outputPath, std::ios::binary);
  if (!output) throw std::runtime_error("cannot create " + outputPath);
  writeValue(output, kMagic);
  writeValue(output, kVersion);
  writeValue(output, nrow);
  writeValue(output, ncol);
  writeValue(output, static_cast<std::int32_t>(meta.nLayer));
  writeValue(output, static_cast<std::int32_t>(meta.nSector));
  writeValue(output, static_cast<std::int32_t>(meta.segmentationMode));
  writeValue(output, static_cast<std::int32_t>(meta.physicsFlag));
  writeValue(output, static_cast<std::int32_t>(meta.nSegTH));
  writeValue(output, static_cast<std::int32_t>(meta.nSegTLC));
  writeValue(output, meta.thetaMinDeg);
  writeValue(output, meta.thetaMaxDeg);
  writeValue(output, threshold);
  for (Long64_t entry = 0; entry < nrow64; ++entry) {
    calarr->GetEntry(entry);
    if (!energy) throw std::runtime_error("null calarr energy vector");
    const auto features = makeFeatures(meta, eventID, *energy, threshold);
    output.write(reinterpret_cast<const char *>(features.data()),
                 static_cast<std::streamsize>(features.size() * sizeof(float)));
  }
  output.close();
  writeManifest(outputPath, inputPath, meta, threshold, nrow, names);
  std::cout << inputPath << " -> " << outputPath << " rows=" << nrow
            << " cols=" << ncol << " geometry=" << meta.nLayer << "x"
            << meta.nSector << " threshold=" << threshold << " MeV\n";
}
} // namespace

int main(int argc, char **argv) {
  if (argc != 4) {
    std::cerr << "usage: bgo_extract_pattern_v1 input.root output.bgo2 threshold_MeV\n";
    return 2;
  }
  try {
    extract(argv[1], argv[2], std::stod(argv[3]));
  } catch (const std::exception &error) {
    std::cerr << "ERROR: " << error.what() << "\n";
    return 1;
  }
  return 0;
}
