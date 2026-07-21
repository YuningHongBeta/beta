#ifndef BGOeggGeometry_h
#define BGOeggGeometry_h 1

#include <cmath>
#include <stdexcept>
#include <vector>

namespace BGOeggGeometry
{
// Ideal no-gap crystal geometry from ELPH Annual Report 2010, pp. 62--65,
// equations (1), (2), (8), and (10).  The 31-ring shape keeps the published
// five forward and four backward extensions.  beta_org is the z-mirror of the
// K18 installation: its +z beam entrance is the long bottom/small opening and
// its -z downstream side is the short head/large opening.
struct Ring
{
  double thetaLowDeg;
  double thetaHighDeg;
  double frontRadiusMm;
  bool forward;
};

constexpr int kSectorCount = 60;
constexpr int kNominalForwardRings = 13;
constexpr int kNominalBackwardRings = 9;
constexpr double kPhiSpanDeg = 6.0;
constexpr double kCrystalLengthMm = 220.0;
constexpr double kForwardDivisionB = 1.3;

inline double Pi()
{
  return std::acos(-1.0);
}

inline double ForwardThetaLowDeg(int type)
{
  if (type < 0)
    throw std::runtime_error("BGOegg forward crystal type must be non-negative");
  const double delta = kPhiSpanDeg * Pi() / 180.0;
  const double theta = Pi() / 2.0 - std::atan(
      kForwardDivisionB * std::tan((type + 1) * delta / kForwardDivisionB));
  return theta * 180.0 / Pi();
}

inline double ForwardThetaHighDeg(int type)
{
  if (type < 0)
    throw std::runtime_error("BGOegg forward crystal type must be non-negative");
  if (type == 0)
    return 90.0;
  return ForwardThetaLowDeg(type - 1);
}

inline double ForwardRadiusMm(double thetaLowDeg, double thetaHighDeg)
{
  const double theta = 0.5 * (thetaLowDeg + thetaHighDeg) * Pi() / 180.0;
  const double denominator = std::sqrt(
      0.25 * std::cos(theta) * std::cos(theta) +
      std::sin(theta) * std::sin(theta));
  return 200.0 / denominator;
}

inline std::vector<Ring> BuildRings(int nLayer)
{
  if (nLayer != 22 && nLayer != 31)
    throw std::runtime_error("BGOegg published geometry supports 22 or 31 layers");

  const int forwardExtra = nLayer == 31 ? 5 : 0;
  const int backwardExtra = nLayer == 31 ? 4 : 0;
  const int nForward = kNominalForwardRings + forwardExtra;
  const int nBackward = kNominalBackwardRings + backwardExtra;

  std::vector<Ring> rings;
  rings.reserve(nForward + nBackward);
  for (int type = nForward - 1; type >= 0; --type)
  {
    const double thetaLow = ForwardThetaLowDeg(type);
    const double thetaHigh = ForwardThetaHighDeg(type);
    rings.push_back(
        {thetaLow, thetaHigh, ForwardRadiusMm(thetaLow, thetaHigh), true});
  }
  for (int type = 0; type < nBackward; ++type)
  {
    const double thetaLow = 90.0 + type * kPhiSpanDeg;
    rings.push_back(
        {thetaLow, thetaLow + kPhiSpanDeg, 200.0, false});
  }
  return rings;
}
} // namespace BGOeggGeometry

#endif
