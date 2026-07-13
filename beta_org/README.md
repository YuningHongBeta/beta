# beta_org

Current Geant4 (v11) simulation for the Lambda beta-decay feasibility study.
Upgraded from Cal_src (v10) with configurable π⁻ absorption physics and more realistic simulation.

## Physics

- **Primary**: Lambda 3-body beta decay (Λ → p e⁻ ν̄) — electron energy spectrum sampled via Dalitz-plot rejection (Källén function)
- **Alternative primaries**: π⁻ or π⁰ at fixed momentum (`/gen/mode e|pim|pi0`)
- **Physics**: tuned INCLXX+GEM baseline (`PhysicsFlag=4`): BGO π⁻
  inelastic XS ×1.65, π⁻-origin neutron yield ×3 at/above 20 MeV and ×0.5
  below 20 MeV, and nuclear-like gamma yield ×0.5. These parameters are fixed
  during geometry/readout scans.
- **Cross-section biasing**: configurable neutron/inelastic bias factors for enhanced statistics

## Detector geometry

- Spherical BGO calorimeter: R_min=30 cm, thickness=20 cm, 15 layers × 15 sectors
- Beam-hole openings: upstream 5.7° (30.7 msr), downstream 9.7° (89.8 msr)
- Enriched ⁶Li target (90%, cylindrical, r=1.5 cm, L=30 cm) at center
- TH hodoscope: plastic scintillator, r=1.5–2.1 cm, 30 segments
- TLC hodoscope: acrylic Cherenkov, r=2.2–2.8 cm, 30 segments

TH and TLC are sensitive detectors. The `th` tree stores per-segment deposited
energy, earliest deposit time, and charged path. The `tlc` tree stores an
analytic Frank–Tamm response for Lucite (`n=1.49`, 300–600 nm): earliest
Cherenkov time, above-threshold path, and produced expected photons. Optical
transport, collection efficiency, and photon-detection efficiency are not yet
applied. TLC deposited energy and total charged path are saved only as
`*_truth` diagnostics and must not be classifier inputs.

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
mkdir -p output

# Reproducible 100k electron sample, current 15x15 uniform-theta geometry
BETA_PRIMARY=e BETA_OUTPUT=output/incl_e_15x15 \
  BETA_SEED=6302026 BETA_THREADS=8 ./beta macros/run_e_100k.mac

# Matching 100k pi-minus sample
BETA_PRIMARY=pim BETA_OUTPUT=output/incl_pim_15x15 \
  BETA_SEED=6302026 BETA_THREADS=8 ./beta macros/run_pim_100k.mac

./beta  # interactive mode with visualization
```

Geometry/readout scan controls are environment variables:

- `BETA_N_LAYER`, `BETA_N_SECTOR`: BGO theta/phi segmentation (default 15x15).
- `BETA_SEGMENTATION`: `uniform_theta` or `equal_solid_angle`.
- `BETA_PRIMARY`: `e`, `pim`, or `pi0`.
- `BETA_OUTPUT`: output path without `.root`.
- `BETA_WRITE_CALHIT=0`: suppress the large per-hit tree during geometry scans;
  the event-level `calarr` vector remains available.
- `BETA_SEED`, `BETA_THREADS`: reproducibility and worker count.

Every output contains a `runmeta` tree with the geometry, primary, seed, and
fixed physics parameters. The tuned INCL physics parameters are not scan axes.
