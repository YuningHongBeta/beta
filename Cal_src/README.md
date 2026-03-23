# Cal_src

Original (Geant4 v10) simulation for the Lambda beta-decay feasibility study. Simulates particles (π⁻ 100 MeV/c, Lambda beta decay electrons) in a spherical BGO calorimeter + ⁶Li target geometry, recording per-crystal energy deposits and hit patterns. Superseded by [beta_org/](../beta_org/) which uses Geant4 v11 with configurable physics.

## Build

Requires Geant4 and ROOT.

```bash
source env_sim.sh
mkdir -p build && cd build
cmake ..
make -j$(nproc)
```

## Run

```bash
cd build
./cal run1.mac      # batch mode
./cal               # interactive mode with visualization
```

## Structure

```
Cal_src/
├── cal.cc                  # Main program
├── CMakeLists.txt          # CMake build configuration
├── include/                # Header files
│   ├── calDetectorConstruction.hh
│   ├── calPrimaryGeneratorAction.hh
│   ├── CalorimeterSD.hh / CalorimeterHit.hh
│   ├── ScintillatorSD.hh / TargetSD.hh
│   └── ...
├── src/                    # Source files
├── *.mac                   # Geant4 macro files
├── plotHisto.C             # ROOT analysis macro
├── plotNtuple.C            # ROOT ntuple plotting
└── scr_checkHitpatPoly.C   # Hit pattern check script
```
