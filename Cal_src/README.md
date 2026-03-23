# Cal_src

Geant4-based calorimeter simulation for the Hyperball-X Ge detector calibration source study.

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
