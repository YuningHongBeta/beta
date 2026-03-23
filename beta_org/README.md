# beta_org

Current Geant4 (v11) simulation for the Lambda beta-decay feasibility study.
Upgraded from Cal_src (v10) with configurable π⁻ absorption physics and more realistic simulation.

## Physics

- **Primary**: Lambda 3-body beta decay (Λ → p e⁻ ν̄) — electron energy spectrum sampled via Dalitz-plot rejection (Källén function)
- **Alternative primaries**: π⁻ or π⁰ at fixed momentum (`/gen/mode e|pim|pi0`)
- **Physics lists**: FTFP_BERT (default), INCLXX+ABLA++, INCLXX+GEM — selectable via `PhysicsFlag`
- **Cross-section biasing**: configurable neutron/inelastic bias factors for enhanced statistics

## Detector geometry

- Spherical BGO calorimeter: R_min=30 cm, thickness=20 cm, 15 layers × 15 sectors
- Beam-hole openings: upstream 5.7° (30.7 msr), downstream 9.7° (89.8 msr)
- Enriched ⁶Li target (90%, cylindrical, r=1.5 cm, L=30 cm) at center
- TH hodoscope: plastic scintillator, r=1.5–2.1 cm, 30 segments
- TLC hodoscope: acrylic Cherenkov, r=2.2–2.8 cm, 30 segments

## Build

Requires Geant4 v11.

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
