#ifndef BetaConfig_h
#define BetaConfig_h 1

#include "BGOeggGeometry.hh"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <string>

class BetaConfig
{
public:
  static const BetaConfig &Instance()
  {
    static const BetaConfig config;
    return config;
  }

  int NLayer() const { return fNLayer; }
  int NSector() const { return fNSector; }
  int NCells() const { return fNLayer * fNSector; }
  double PhiWidthDeg() const { return 360.0 / static_cast<double>(fNSector); }
  bool EqualSolidAngle() const { return fSegmentation == "equal_solid_angle"; }
  bool BGOeggPublishedSegmentation() const { return fSegmentation == "bgoegg_published"; }
  int SegmentationMode() const
  {
    if (BGOeggPublishedSegmentation())
      return 2;
    return EqualSolidAngle() ? 1 : 0;
  }
  const std::string &Segmentation() const { return fSegmentation; }
  const std::string &Output() const { return fOutput; }
  const std::string &Primary() const { return fPrimary; }
  bool BeamOverlay() const
  {
    return fPrimary == "beam" ||
           fPrimary == "e_beam" ||
           fPrimary == "pim_beam" ||
           fPrimary == "pi0_beam";
  }
  bool BeamOnly() const { return fPrimary == "beam"; }
  const std::string &BeamParticle() const { return fBeamParticle; }
  double BeamMomentumMeVC() const { return fBeamMomentumMeVC; }
  const std::string &BeamMultiplicityMode() const { return fBeamMultiplicityMode; }
  int BeamFixedMultiplicity() const { return fBeamFixedMultiplicity; }
  double BeamMeanPerBgoGate() const { return fBeamMeanPerBgoGate; }
  double BgoGateWidthNs() const { return fBgoGateWidthNs; }
  double HodoGateWidthNs() const { return fHodoGateWidthNs; }
  const std::string &BeamProfileModel() const { return fBeamProfileModel; }
  double BeamXMeanMm() const { return fBeamXMeanMm; }
  double BeamYMeanMm() const { return fBeamYMeanMm; }
  double BeamXSigmaMm() const { return fBeamXSigmaMm; }
  double BeamYSigmaMm() const { return fBeamYSigmaMm; }
  double BeamXMaxAbsMm() const { return fBeamXMaxAbsMm; }
  double BeamYMaxAbsMm() const { return fBeamYMaxAbsMm; }
  double PiMinusMomentumMeVC() const { return fPiMinusMomentumMeVC; }
  double PiZeroMomentumMeVC() const { return fPiZeroMomentumMeVC; }
  const std::string &TargetMaterial() const { return fTargetMaterial; }
  double TargetArealDensityGCM2() const { return fTargetArealDensityGCM2; }
  double TargetDensityGCM3() const { return fTargetDensityGCM3; }
  double TargetRadiusMm() const { return fTargetRadiusMm; }
  double TargetLengthMm() const
  {
    return 10.0 * fTargetArealDensityGCM2 / fTargetDensityGCM3;
  }
  bool WriteCalHit() const { return fWriteCalHit; }
  int Threads() const { return fThreads; }
  long Seed() const { return fSeed; }
  const std::string &Geometry() const { return fGeometry; }
  int GeometryMode() const
  {
    if (fGeometry == "bgoegg_envelope")
      return 1;
    return fGeometry == "bgoegg_frustum" ? 2 : 0;
  }
  bool BGOeggEnvelope() const { return GeometryMode() == 1; }
  bool BGOeggFrustum() const { return GeometryMode() == 2; }
  bool BGOeggGeometry() const { return GeometryMode() != 0; }
  int BGOeggForwardExtraLayers() const { return fNLayer == 31 ? 5 : 0; }
  int BGOeggBackwardExtraLayers() const { return fNLayer == 31 ? 4 : 0; }
  const std::string &PhotonCounter() const { return fPhotonCounter; }
  int PhotonCounterMode() const
  {
    if (fPhotonCounter == "downstream")
      return 1;
    if (fPhotonCounter == "two_sided")
      return 2;
    return fPhotonCounter == "upstream" ? 3 : 0;
  }
  bool HasDownstreamPhotonCounter() const
  {
    return PhotonCounterMode() == 1 || PhotonCounterMode() == 2;
  }
  bool HasUpstreamPhotonCounter() const
  {
    return PhotonCounterMode() == 2 || PhotonCounterMode() == 3;
  }
  int PCNLayers() const { return fPCNLayers; }
  double PCPbThicknessMm() const { return fPCPbThicknessMm; }
  double PCScintiThicknessMm() const { return fPCScintiThicknessMm; }
  double PCZFrontCm() const { return fPCZFrontCm; }
  double PCDownThetaInnerDeg() const { return fPCDownThetaInnerDeg; }
  double PCDownThetaOuterDeg() const { return fPCDownThetaOuterDeg; }
  double PCUpThetaInnerDeg() const { return fPCUpThetaInnerDeg; }
  double PCUpThetaOuterDeg() const { return fPCUpThetaOuterDeg; }
  bool BgoZOffsetConfigured() const { return fBgoZOffsetConfigured; }
  double BgoZOffsetCm() const { return fBgoZOffsetCm; }
  double RMinCm() const { return BGOeggGeometry() ? 20.0 : 30.0; }
  double ThicknessCm() const { return BGOeggGeometry() ? 22.0 : 20.0; }
  double ThetaMinDeg() const
  {
    if (!BGOeggFrustum())
      return fThetaMinDeg;
    return BGOeggGeometry::BuildRings(fNLayer).front().thetaLowDeg;
  }
  double ThetaMaxDeg() const
  {
    if (!BGOeggFrustum())
      return fThetaMaxDeg;
    return BGOeggGeometry::BuildRings(fNLayer).back().thetaHighDeg;
  }
  const char *GeometryModel() const
  {
    if (BGOeggEnvelope())
      return "spherical_shell_envelope_not_exact_bgoegg_trapezoids";
    if (BGOeggFrustum())
      return "published_bgoegg_frusta_ideal_no_gaps";
    return "spherical_shell_current";
  }

private:
  static int ReadInt(const char *name, int fallback, int minimum, int maximum)
  {
    const char *raw = std::getenv(name);
    if (!raw || !*raw)
      return fallback;
    char *end = nullptr;
    const long value = std::strtol(raw, &end, 10);
    if (!end || *end != 0 || value < minimum || value > maximum)
      throw std::runtime_error(std::string("Invalid ") + name + "=" + raw);
    return static_cast<int>(value);
  }

  static double ReadDouble(const char *name, double fallback,
                           double minimum, double maximum)
  {
    const char *raw = std::getenv(name);
    if (!raw || !*raw)
      return fallback;
    char *end = nullptr;
    const double value = std::strtod(raw, &end);
    if (!end || *end != 0 || !std::isfinite(value) ||
        value < minimum || value > maximum)
      throw std::runtime_error(std::string("Invalid ") + name + "=" + raw);
    return value;
  }

  static bool ReadBool(const char *name, bool fallback)
  {
    const char *raw = std::getenv(name);
    if (!raw || !*raw)
      return fallback;
    const std::string value(raw);
    if (value == "1" || value == "true" || value == "TRUE" || value == "yes")
      return true;
    if (value == "0" || value == "false" || value == "FALSE" || value == "no")
      return false;
    throw std::runtime_error(std::string("Invalid ") + name + "=" + raw);
  }

  static std::string ReadString(const char *name, const char *fallback)
  {
    const char *raw = std::getenv(name);
    return (raw && *raw) ? std::string(raw) : std::string(fallback);
  }

  BetaConfig()
      : fGeometry(ReadString("BETA_GEOMETRY", "current")),
        fPhotonCounter(ReadString("BETA_PHOTON_COUNTER", "none")),
        fPCNLayers(ReadInt("BETA_PC_N_LAYERS", 8, 1, 64)),
        fPCPbThicknessMm(ReadDouble("BETA_PC_PB_THICKNESS_MM", 1.0, 0.01, 20.0)),
        fPCScintiThicknessMm(ReadDouble(
            "BETA_PC_SCINTI_THICKNESS_MM", 5.0, 0.1, 100.0)),
        fPCZFrontCm(ReadDouble("BETA_PC_Z_FRONT_CM", 52.0, 1.0, 89.0)),
        fPCDownThetaInnerDeg(ReadDouble(
            "BETA_PC_DOWN_THETA_INNER_DEG", 9.698, 0.01, 80.0)),
        fPCDownThetaOuterDeg(ReadDouble(
            "BETA_PC_DOWN_THETA_OUTER_DEG", 24.0, 0.02, 85.0)),
        fPCUpThetaInnerDeg(ReadDouble(
            "BETA_PC_UP_THETA_INNER_DEG", 5.666, 0.01, 80.0)),
        fPCUpThetaOuterDeg(ReadDouble(
            "BETA_PC_UP_THETA_OUTER_DEG", 36.0, 0.02, 85.0)),
        fBgoZOffsetConfigured(std::getenv("BETA_BGO_Z_OFFSET_CM") != nullptr),
        fBgoZOffsetCm(ReadDouble("BETA_BGO_Z_OFFSET_CM", 0.0, -10.0, 10.0)),
        fThetaMinConfigured(std::getenv("BETA_BGO_THETA_MIN_DEG") != nullptr),
        fThetaMaxConfigured(std::getenv("BETA_BGO_THETA_MAX_DEG") != nullptr),
        fThetaMinDeg(ReadDouble("BETA_BGO_THETA_MIN_DEG",
                               fGeometry == "bgoegg_envelope" ? 24.0 : 5.666,
                               0.1, 179.8)),
        fThetaMaxDeg(ReadDouble("BETA_BGO_THETA_MAX_DEG",
                               fGeometry == "bgoegg_envelope" ? 144.0 : 170.302,
                               0.2, 179.9)),
        fNLayer(ReadInt("BETA_N_LAYER", (fGeometry == "bgoegg_envelope" ||
                                         fGeometry == "bgoegg_frustum") ? 22 : 15,
                        1, 200)),
        fNSector(ReadInt("BETA_N_SECTOR", (fGeometry == "bgoegg_envelope" ||
                                          fGeometry == "bgoegg_frustum") ? 60 : 15,
                         1, 360)),
        fSegmentation(ReadString(
            "BETA_SEGMENTATION",
            fGeometry == "bgoegg_frustum" ? "bgoegg_published" : "uniform_theta")),
        fOutput(ReadString("BETA_OUTPUT", "output/beta")),
        fPrimary(ReadString("BETA_PRIMARY", "e")),
        fBeamParticle(ReadString("BETA_BEAM_PARTICLE", "pi+")),
        fBeamMomentumMeVC(ReadDouble(
            "BETA_BEAM_MOMENTUM_MEV_C", 1100.0, 1.0, 10000.0)),
        fBeamMultiplicityMode(ReadString("BETA_BEAM_MULTIPLICITY_MODE", "fixed")),
        fBeamFixedMultiplicity(ReadInt("BETA_BEAM_FIXED_MULTIPLICITY", 1, 0, 100000)),
        fBeamMeanPerBgoGate(ReadDouble(
            "BETA_BEAM_MEAN_PER_BGO_GATE", 1.0, 0.0, 100000.0)),
        fBgoGateWidthNs(ReadDouble("BETA_BGO_GATE_WIDTH_NS", 0.0, 0.0, 1.0e9)),
        fHodoGateWidthNs(ReadDouble("BETA_HODO_GATE_WIDTH_NS", 0.0, 0.0, 1.0e9)),
        fBeamProfileModel(ReadString("BETA_BEAM_PROFILE_MODEL", "pencil")),
        fBeamXMeanMm(ReadDouble("BETA_BEAM_X_MEAN_MM", 0.0, -1000.0, 1000.0)),
        fBeamYMeanMm(ReadDouble("BETA_BEAM_Y_MEAN_MM", 0.0, -1000.0, 1000.0)),
        fBeamXSigmaMm(ReadDouble("BETA_BEAM_X_SIGMA_MM", 0.0, 0.0, 1000.0)),
        fBeamYSigmaMm(ReadDouble("BETA_BEAM_Y_SIGMA_MM", 0.0, 0.0, 1000.0)),
        fBeamXMaxAbsMm(ReadDouble("BETA_BEAM_X_MAX_ABS_MM", 0.0, 0.0, 1000.0)),
        fBeamYMaxAbsMm(ReadDouble("BETA_BEAM_Y_MAX_ABS_MM", 0.0, 0.0, 1000.0)),
        fPiMinusMomentumMeVC(ReadDouble(
            "BETA_PIM_MOMENTUM_MEV_C", 100.0, 0.01, 1000.0)),
        fPiZeroMomentumMeVC(ReadDouble(
            "BETA_PI0_MOMENTUM_MEV_C", 100.0, 0.01, 1000.0)),
        fTargetMaterial(ReadString("BETA_TARGET_MATERIAL", "Li6_90pct")),
        fTargetArealDensityGCM2(ReadDouble(
            "BETA_TARGET_AREAL_DENSITY_G_CM2", 14.1, 0.001, 100.0)),
        fTargetDensityGCM3(ReadDouble(
            "BETA_TARGET_DENSITY_G_CM3", 0.47, 0.001, 30.0)),
        fTargetRadiusMm(ReadDouble(
            "BETA_TARGET_RADIUS_MM", 15.0, 0.1, 15.0)),
        fWriteCalHit(ReadBool("BETA_WRITE_CALHIT", true)),
        fThreads(ReadInt("BETA_THREADS", 8, 1, 256)),
        fSeed(ReadInt("BETA_SEED", 6302026, 1, 2147483647))
  {
    if (fGeometry != "current" && fGeometry != "bgoegg_envelope" &&
        fGeometry != "bgoegg_frustum")
      throw std::runtime_error(
          "BETA_GEOMETRY must be current, bgoegg_envelope, or bgoegg_frustum");
    if (fPhotonCounter != "none" && fPhotonCounter != "downstream" &&
        fPhotonCounter != "upstream" && fPhotonCounter != "two_sided")
      throw std::runtime_error(
          "BETA_PHOTON_COUNTER must be none, downstream, upstream, or two_sided");
    if (fPCPbThicknessMm == 0.0 && fPCScintiThicknessMm == 0.0)
      throw std::runtime_error("Photon counter requires non-zero layer thickness");
    if (fPCDownThetaInnerDeg >= fPCDownThetaOuterDeg ||
        fPCUpThetaInnerDeg >= fPCUpThetaOuterDeg)
      throw std::runtime_error(
          "Photon-counter inner angles must be smaller than outer angles");
    const double pcBackCm = fPCZFrontCm +
        0.1 * fPCNLayers * (fPCPbThicknessMm + fPCScintiThicknessMm);
    if (pcBackCm >= 90.0)
      throw std::runtime_error("Photon counter exceeds the world z boundary");
    const double pcMaxRadiusCm = pcBackCm * std::tan(
        std::max(fPCDownThetaOuterDeg, fPCUpThetaOuterDeg) *
        BGOeggGeometry::Pi() / 180.0);
    if (pcMaxRadiusCm >= 90.0)
      throw std::runtime_error("Photon counter exceeds the world radial boundary");
    if (BGOeggFrustum() && fPhotonCounter != "none")
    {
      double downstreamMaxCm = -1.0e9;
      double upstreamExtentCm = -1.0e9;
      for (const auto &ring : BGOeggGeometry::BuildRings(fNLayer))
      {
        for (const double radiusMm : {
                 ring.frontRadiusMm,
                 ring.frontRadiusMm + BGOeggGeometry::kCrystalLengthMm})
        {
          for (const double thetaDeg : {ring.thetaLowDeg, ring.thetaHighDeg})
          {
            const double zCm = fBgoZOffsetCm + 0.1 * radiusMm * std::cos(
                thetaDeg * BGOeggGeometry::Pi() / 180.0);
            downstreamMaxCm = std::max(downstreamMaxCm, zCm);
            upstreamExtentCm = std::max(upstreamExtentCm, -zCm);
          }
        }
      }
      constexpr double clearanceCm = 0.05;
      if (HasDownstreamPhotonCounter() &&
          fPCZFrontCm <= downstreamMaxCm + clearanceCm)
        throw std::runtime_error(
            "Downstream photon counter intersects the published BGOegg frusta");
      if (HasUpstreamPhotonCounter() &&
          fPCZFrontCm <= upstreamExtentCm + clearanceCm)
        throw std::runtime_error(
            "Upstream photon counter intersects the published BGOegg frusta");
    }
    if (!BGOeggGeometry() && fPhotonCounter != "none" &&
        fPCZFrontCm <= RMinCm() + ThicknessCm() + 0.05)
      throw std::runtime_error(
          "Photon counter intersects the current spherical BGOC");
    if (!BGOeggGeometry() && fBgoZOffsetCm != 0.0)
      throw std::runtime_error(
          "BETA_BGO_Z_OFFSET_CM requires a BGOegg geometry");
    if (BGOeggFrustum() &&
        (fThetaMinConfigured || fThetaMaxConfigured))
      throw std::runtime_error(
          "BGO theta overrides are not supported for published BGOegg frusta");
    if (fThetaMinDeg >= fThetaMaxDeg)
      throw std::runtime_error(
          "BETA_BGO_THETA_MIN_DEG must be less than BETA_BGO_THETA_MAX_DEG");
    if (fSegmentation != "uniform_theta" &&
        fSegmentation != "equal_solid_angle" &&
        fSegmentation != "bgoegg_published")
      throw std::runtime_error(
          "BETA_SEGMENTATION must be uniform_theta, equal_solid_angle, or bgoegg_published");
    if (BGOeggFrustum() && fSegmentation != "bgoegg_published")
      throw std::runtime_error(
          "BETA_GEOMETRY=bgoegg_frustum requires BETA_SEGMENTATION=bgoegg_published");
    if (!BGOeggFrustum() && fSegmentation == "bgoegg_published")
      throw std::runtime_error(
          "BETA_SEGMENTATION=bgoegg_published requires BETA_GEOMETRY=bgoegg_frustum");
    if (BGOeggFrustum() && fNSector != 60)
      throw std::runtime_error("BGOegg published frusta require BETA_N_SECTOR=60");
    if (BGOeggFrustum() && fNLayer != 22 && fNLayer != 31)
      throw std::runtime_error("BGOegg published frusta support BETA_N_LAYER=22 or 31");
    if (fPrimary != "e" && fPrimary != "pim" && fPrimary != "pi0" &&
        fPrimary != "beam" && fPrimary != "e_beam" &&
        fPrimary != "pim_beam" && fPrimary != "pi0_beam")
      throw std::runtime_error(
          "BETA_PRIMARY must be e, pim, pi0, beam, e_beam, pim_beam, or pi0_beam");
    if (fBeamParticle != "pi+" && fBeamParticle != "pi-" &&
        fBeamParticle != "kaon+" && fBeamParticle != "kaon-")
      throw std::runtime_error(
          "BETA_BEAM_PARTICLE must be pi+, pi-, kaon+, or kaon-");
    if (fBeamMultiplicityMode != "fixed" && fBeamMultiplicityMode != "poisson")
      throw std::runtime_error(
          "BETA_BEAM_MULTIPLICITY_MODE must be fixed or poisson");
    if (fBeamMultiplicityMode == "poisson" && fBgoGateWidthNs <= 0.0)
      throw std::runtime_error(
          "BETA_BGO_GATE_WIDTH_NS must be positive for poisson beam multiplicity");
    if (fBeamProfileModel != "pencil" &&
        fBeamProfileModel != "independent_truncated_gaussian")
      throw std::runtime_error(
          "BETA_BEAM_PROFILE_MODEL must be pencil or independent_truncated_gaussian");
    if (fBeamProfileModel == "independent_truncated_gaussian" &&
        (fBeamXSigmaMm <= 0.0 || fBeamYSigmaMm <= 0.0 ||
         fBeamXMaxAbsMm <= 0.0 || fBeamYMaxAbsMm <= 0.0 ||
         std::abs(fBeamXMeanMm) >= fBeamXMaxAbsMm ||
         std::abs(fBeamYMeanMm) >= fBeamYMaxAbsMm))
      throw std::runtime_error(
          "truncated Gaussian beam requires positive sigma/max and |mean| < max");
    if (fTargetMaterial != "Li6_90pct" &&
        fTargetMaterial != "C13_100pct")
      throw std::runtime_error(
          "BETA_TARGET_MATERIAL must be Li6_90pct or C13_100pct");
    if (TargetLengthMm() >= 600.0)
      throw std::runtime_error(
          "Configured target length must be smaller than the 600 mm TH/TLC length");
    if (NCells() > 20000)
      throw std::runtime_error("BETA_N_LAYER * BETA_N_SECTOR exceeds 20000");
  }

  std::string fGeometry;
  std::string fPhotonCounter;
  int fPCNLayers;
  double fPCPbThicknessMm;
  double fPCScintiThicknessMm;
  double fPCZFrontCm;
  double fPCDownThetaInnerDeg;
  double fPCDownThetaOuterDeg;
  double fPCUpThetaInnerDeg;
  double fPCUpThetaOuterDeg;
  bool fBgoZOffsetConfigured;
  double fBgoZOffsetCm;
  bool fThetaMinConfigured;
  bool fThetaMaxConfigured;
  double fThetaMinDeg;
  double fThetaMaxDeg;
  int fNLayer;
  int fNSector;
  std::string fSegmentation;
  std::string fOutput;
  std::string fPrimary;
  std::string fBeamParticle;
  double fBeamMomentumMeVC;
  std::string fBeamMultiplicityMode;
  int fBeamFixedMultiplicity;
  double fBeamMeanPerBgoGate;
  double fBgoGateWidthNs;
  double fHodoGateWidthNs;
  std::string fBeamProfileModel;
  double fBeamXMeanMm;
  double fBeamYMeanMm;
  double fBeamXSigmaMm;
  double fBeamYSigmaMm;
  double fBeamXMaxAbsMm;
  double fBeamYMaxAbsMm;
  double fPiMinusMomentumMeVC;
  double fPiZeroMomentumMeVC;
  std::string fTargetMaterial;
  double fTargetArealDensityGCM2;
  double fTargetDensityGCM3;
  double fTargetRadiusMm;
  bool fWriteCalHit;
  int fThreads;
  long fSeed;
};

#endif
