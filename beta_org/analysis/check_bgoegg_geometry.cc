#include "../include/BGOeggGeometry.hh"

#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
void requireClose(double actual, double expected, double tolerance,
                  const char *label)
{
  if (std::abs(actual - expected) > tolerance)
    throw std::runtime_error(
        std::string(label) + ": actual=" + std::to_string(actual) +
        " expected=" + std::to_string(expected));
}
} // namespace

int main()
{
  try
  {
    const auto nominal = BGOeggGeometry::BuildRings(22);
    const auto extended = BGOeggGeometry::BuildRings(31);
    if (nominal.size() != 22 || extended.size() != 31)
      throw std::runtime_error("unexpected BGOegg ring count");

    constexpr std::array<double, 14> forwardEdges = {
        23.95, 27.97, 32.13, 36.46, 40.97, 45.67, 50.58,
        55.69, 61.02, 66.54, 72.23, 78.07, 84.01, 90.00};
    constexpr std::array<double, 13> forwardRadii = {
        318.76, 302.18, 286.27, 271.32, 257.56, 245.11, 234.07,
        224.51, 216.47, 210.00, 205.11, 201.84, 200.21};
    for (std::size_t ring = 0; ring < forwardRadii.size(); ++ring)
    {
      requireClose(nominal[ring].thetaLowDeg, forwardEdges[ring], 0.015,
                   "nominal forward thetaLow");
      requireClose(nominal[ring].thetaHighDeg, forwardEdges[ring + 1], 0.015,
                   "nominal forward thetaHigh");
      requireClose(nominal[ring].frontRadiusMm, forwardRadii[ring], 0.015,
                   "nominal forward radius");
    }
    for (int ring = 0; ring < 9; ++ring)
    {
      const auto &entry = nominal[13 + ring];
      requireClose(entry.thetaLowDeg, 90.0 + 6.0 * ring, 1.0e-12,
                   "nominal backward thetaLow");
      requireClose(entry.thetaHighDeg, 96.0 + 6.0 * ring, 1.0e-12,
                   "nominal backward thetaHigh");
      requireClose(entry.frontRadiusMm, 200.0, 1.0e-12,
                   "nominal backward radius");
    }

    requireClose(extended.front().thetaLowDeg, 5.336032242257286, 1.0e-12,
                 "extended thetaMin");
    requireClose(extended.back().thetaHighDeg, 168.0, 1.0e-12,
                 "extended thetaMax");
    const double coverage = 0.5 * (
        std::cos(extended.front().thetaLowDeg * BGOeggGeometry::Pi() / 180.0) -
        std::cos(extended.back().thetaHighDeg * BGOeggGeometry::Pi() / 180.0));

    std::cout << std::fixed << std::setprecision(6)
              << "BGOegg geometry check OK\n"
              << "nominal: rings=22 cells=1320 theta="
              << nominal.front().thetaLowDeg << "--"
              << nominal.back().thetaHighDeg << " deg\n"
              << "extended: rings=31 cells=1860 theta="
              << extended.front().thetaLowDeg << "--"
              << extended.back().thetaHighDeg << " deg coverage="
              << 100.0 * coverage << "% of 4pi\n";
    return 0;
  }
  catch (const std::exception &error)
  {
    std::cerr << "ERROR: " << error.what() << '\n';
    return 1;
  }
}
