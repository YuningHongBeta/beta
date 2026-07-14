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
- TH hodoscope: plastic scintillator, r=1.5–2.1 cm, z=-300–+300 mm,
  30 segments
- TLC hodoscope: acrylic Cherenkov, r=2.2–2.8 cm, 30 segments

For feasibility scans, `BETA_GEOMETRY=bgoegg_envelope` selects a spherical
shell envelope with R_min=20 cm, thickness=22 cm, theta=24–144 degrees, and
default segmentation 22x60. It reproduces the published active depth and polar
coverage only. It does **not** model the exact BGOegg trapezoidal crystals,
forward prolate-spheroid inner face, support, gaps, PMTs, or services, and must
not be used as an engineering geometry.

`BETA_PHOTON_COUNTER=downstream|two_sided` adds provisional unsegmented collars
for acceptance studies. Each collar is `8 x (1 mm Pb + 5 mm plastic)` beginning
at `|z|=52 cm`. The downstream (+z) collar covers 9.698–24 degrees and the
optional upstream (-z) collar covers 5.666–36 degrees. `evt.EdepPC_MeV` is the
sum of energy deposited in all enabled plastic layers. Optical response,
threshold, segmentation, timing, pile-up, supports, and readout are absent.

TH and TLC are sensitive detectors. The `th` tree stores per-segment deposited
energy, earliest step-midpoint deposit time, and analytic arrival times at both
TH ends. The left MPPC is the z=-300 mm end and the right MPPC is the z=+300 mm
end. Its event-level vectors are `dE_MeV`, `time_ns`, `timeLeft_ns`,
`timeRight_ns`, `timeLeftMinusRight_ns`, and `zReco_mm`. With the default
effective light speed `v_eff=150 mm/ns`,

```text
t_left  = t_deposit + (z_deposit - (-300 mm)) / v_eff
t_right = t_deposit + ((+300 mm) - z_deposit) / v_eff
dt_left_minus_right = t_left - t_right
z_reco = 0.5 * v_eff * dt_left_minus_right
```

Thus positive `timeLeftMinusRight_ns` gives positive `zReco_mm`. Each endpoint
time is the earliest analytic arrival over energy-depositing steps in that
segment. The response model is `analytic_step_midpoint_earliest_no_smearing`:
there is no timing offset, time walk, attenuation, threshold, photon
statistics, or smearing. The effective speed is an explicit provisional
detector-response assumption, not a measured TH calibration.
Missing endpoint-derived `timeLeftMinusRight_ns` and `zReco_mm` values use the
out-of-range sentinel `-9999`; endpoint arrival-time sentinels are `-1 ns`.
The `runmeta` tree records the timing assumptions together with the TH shell
radii (`thRMin_mm=15`, `thRMax_mm=21`) needed for geometry-driven path length.

For analysis, TH `dE/dx` and the residual `Delta-z` must use the geometric path
and z intersection of the target-to-leading-BGO-cell line with the TH shell.
`chargedPath_truth_mm` is retained only as a simulation diagnostic and must not
be used as a classifier feature.

The `tlc` tree stores an analytic Frank–Tamm response for Lucite (`n=1.49`,
300–600 nm): earliest
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

- `BETA_GEOMETRY`: `current` (default) or the provisional
  `bgoegg_envelope`. The latter defaults to 22x60 cells.
- `BETA_PHOTON_COUNTER`: `none` (default), `downstream`, or `two_sided`.
- `BETA_BGO_Z_OFFSET_CM`: optional BGOegg-only beam-axis translation in
  `[-10,10]` cm. The target, TH, TLC, and photon counters remain fixed.
- `BETA_N_LAYER`, `BETA_N_SECTOR`: BGO theta/phi segmentation (default 15x15).
- `BETA_SEGMENTATION`: `uniform_theta` or `equal_solid_angle`.
- `BETA_PRIMARY`: `e`, `pim`, or `pi0`.
- `BETA_OUTPUT`: output path without `.root`.
- `BETA_WRITE_CALHIT=0`: suppress the large per-hit tree during geometry scans;
  the event-level `calarr` vector remains available.
- `BETA_SEED`, `BETA_THREADS`: reproducibility and worker count.

Every output contains a `runmeta` tree with the geometry/counter modes and
dimensions, primary, seed, and fixed physics parameters. The tuned INCL physics
parameters are not scan axes.

Example BGOegg-envelope run:

```bash
BETA_GEOMETRY=bgoegg_envelope BETA_PHOTON_COUNTER=downstream \
  BETA_PRIMARY=pi0 BETA_OUTPUT=output/egg_down_pi0 \
  ./beta macros/run.mac
```
BGOegg envelope coverage scans may override the active polar-angle range with
`BETA_BGO_THETA_MIN_DEG` and `BETA_BGO_THETA_MAX_DEG`.  These options are for
the provisional spherical envelope only; they do not model exact BGOegg
trapezoids, supports, services, or readout clearances.
