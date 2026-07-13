#ifndef BetaConfig_h
#define BetaConfig_h 1

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
  int SegmentationMode() const { return EqualSolidAngle() ? 1 : 0; }
  const std::string &Segmentation() const { return fSegmentation; }
  const std::string &Output() const { return fOutput; }
  const std::string &Primary() const { return fPrimary; }
  bool WriteCalHit() const { return fWriteCalHit; }
  int Threads() const { return fThreads; }
  long Seed() const { return fSeed; }

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
      : fNLayer(ReadInt("BETA_N_LAYER", 15, 1, 200)),
        fNSector(ReadInt("BETA_N_SECTOR", 15, 1, 360)),
        fSegmentation(ReadString("BETA_SEGMENTATION", "uniform_theta")),
        fOutput(ReadString("BETA_OUTPUT", "output/beta")),
        fPrimary(ReadString("BETA_PRIMARY", "e")),
        fWriteCalHit(ReadBool("BETA_WRITE_CALHIT", true)),
        fThreads(ReadInt("BETA_THREADS", 8, 1, 256)),
        fSeed(ReadInt("BETA_SEED", 6302026, 1, 2147483647))
  {
    if (fSegmentation != "uniform_theta" &&
        fSegmentation != "equal_solid_angle")
      throw std::runtime_error("BETA_SEGMENTATION must be uniform_theta or equal_solid_angle");
    if (fPrimary != "e" && fPrimary != "pim" && fPrimary != "pi0")
      throw std::runtime_error("BETA_PRIMARY must be e, pim, or pi0");
    if (NCells() > 20000)
      throw std::runtime_error("BETA_N_LAYER * BETA_N_SECTOR exceeds 20000");
  }

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