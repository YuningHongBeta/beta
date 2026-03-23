# beta_org

Geant4 simulation of beta-ray detectors with custom physics list and biasing support.

## Build

Requires Geant4.

```bash
mkdir -p build && cd build
cmake ..
make -j$(nproc)
```

## Run

```bash
cd build
./beta macros/run.mac   # batch mode
./beta                  # interactive mode with visualization
```

## Structure

```
beta_org/
├── src/main.cc                        # Main program
├── CMakeLists.txt                     # CMake build configuration
├── include/                           # Header files
│   ├── betaDetectorConstruction.hh    # Geometry definition
│   ├── betaPhysicsList.hh            # Custom physics list
│   ├── betaPrimaryGeneratorAction.hh  # Particle gun
│   ├── betaBiasingOperator.hh         # Cross-section biasing
│   ├── CalorimeterSD.hh / TargetSD.hh
│   └── ...
├── src/                               # Source files
└── macros/                            # Geant4 macro files
```
