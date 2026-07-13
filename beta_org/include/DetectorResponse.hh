#ifndef DetectorResponse_h
#define DetectorResponse_h 1

#include "G4SystemOfUnits.hh"
#include "globals.hh"

namespace DetectorResponse
{
constexpr G4double THBarZMin = -300.0 * mm;
constexpr G4double THBarZMax = 300.0 * mm;
constexpr G4double THEffectiveLightSpeed = 150.0 * mm / ns;
constexpr G4bool THTimingSmearingApplied = false;
constexpr const char THTimingModel[] =
    "analytic_step_midpoint_earliest_no_smearing";
} // namespace DetectorResponse

#endif