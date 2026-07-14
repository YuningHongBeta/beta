#ifndef BetaConfig_h
#define BetaConfig_h 1

#include "BGOeggGeometry.hh"

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
    return fPhotonCounter == "two_sided" ? 2 : 0;
  }
  bool HasDownstreamPhotonCounter() const { return PhotonCounterMode() >= 1; }
  bool HasUpstreamPhotonCounter() const { return PhotonCounterMode() == 2; }
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
        fWriteCalHit(ReadBool("BETA_WRITE_CALHIT", true)),
        fThreads(ReadInt("BETA_THREADS", 8, 1, 256)),
        fSeed(ReadInt("BETA_SEED", 6302026, 1, 2147483647))
  {
    if (fGeometry != "current" && fGeometry != "bgoegg_envelope" &&
        fGeometry != "bgoegg_frustum")
      throw std::runtime_error(
          "BETA_GEOMETRY must be current, bgoegg_envelope, or bgoegg_frustum");
    if (fPhotonCounter != "none" && fPhotonCounter != "downstream" &&
        fPhotonCounter != "two_sided")
      throw std::runtime_error("BETA_PHOTON_COUNTER must be none, downstream, or two_sided");
    if (!BGOeggGeometry() && fPhotonCounter != "none")
      throw std::runtime_error(
          "BETA_PHOTON_COUNTER collars require a BGOegg geometry");
    if (!BGOeggGeometry() && fBgoZOffsetCm != 0.0)
      throw std::runtime_error(
          "BETA_BGO_Z_OFFSET_CM requires a BGOegg geometry");
    if (!BGOeggEnvelope() &&
        (fThetaMinConfigured || fThetaMaxConfigured))
      throw std::runtime_error(
          "BGO theta overrides require BETA_GEOMETRY=bgoegg_envelope");
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
    if (fPrimary != "e" && fPrimary != "pim" && fPrimary != "pi0")
      throw std::runtime_error("BETA_PRIMARY must be e, pim, or pi0");
    if (NCells() > 20000)
      throw std::runtime_error("BETA_N_LAYER * BETA_N_SECTOR exceeds 20000");
  }

  std::string fGeometry;
  std::string fPhotonCounter;
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
  bool fWriteCalHit;
  int fThreads;
  long fSeed;
};

#endif
