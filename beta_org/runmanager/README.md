# beta_org small runmanager

This directory provides a deliberately small LSF wrapper for reproducible BGO
geometry/readout samples. It calls `scripts/run_bgo_sample.sh`; it does not copy
the E63 analyzer runmanager and it does not change physics.

## Manifest contract

The accepted schemas are enumerated in `runmanager.py`; unknown manifest fields,
unknown schemas, modified manifests with existing state, unsupported geometry,
and event counts outside the checked-in 20k/100k/2M macros are rejected.
`beta-bgo-th-tlc-beam-overlay-v6` adds explicit `beam`, `target`, and `signal`
mappings and accepts clean, beam-only, and same-event overlay primaries.
`beta-bgo-th-tlc-rate-overlay-v7` adds fixed/Poisson beam multiplicity, BGO and
TH/TLC gates, and pencil/truncated-Gaussian target-plane profiles. A list in
`beam.hodo_gate_width_ns` expands a gate scan in addition to the geometry and
primary product. The supplied Kamada K1.1 and Hong K1.8 examples scan
4/10/20/50 ns. These are provisional sensitivity points, not measured beta
electronics gates. Each is a full-width rectangular window centered at the
target signal `t=0`, so time of flight is part of the acceptance.

Relative `build_dir` paths are resolved from `beta_org/`, not from the manifest
directory. `geometries` and `primaries` form a Cartesian product. The current
geometry/timing example is
`examples/current_geometry_timing.yml`. The BGOegg envelope/counter matrix is
`examples/bgoegg_feasibility.yml`.

Each job invokes

```text
scripts/run_bgo_sample.sh RUN_TAG N_LAYER N_SECTOR SEGMENTATION PRIMARY GEOMETRY PHOTON_COUNTER
```

with `BETA_BUILD_DIR`, `BETA_THREADS`, and `BETA_SEED` set from the manifest.
The runner disables the large `calhit` payload and writes
`BUILD_DIR/output/scan/RUN_TAG_PRIMARY.root`.

## Commands

Inspect every `bsub` command without creating or modifying state and without
submitting jobs:

```bash
python3 runmanager/runmanager.py submit \
  runmanager/examples/current_geometry_timing.yml --dry-run
```

This prints the complete manifest expansion even when state already exists.
Use `--retry --dry-run` instead to inspect only currently eligible
failed/missing retry commands.

Submit the full geometry x primary product:

```bash
python3 runmanager/runmanager.py submit \
  runmanager/examples/current_geometry_timing.yml
```

Query LSF state. This command only calls `bjobs`; it never calls `bkill`:

```bash
python3 runmanager/runmanager.py status \
  runmanager/examples/current_geometry_timing.yml
```

Validate completed outputs:

```bash
python3 runmanager/runmanager.py validate \
  runmanager/examples/current_geometry_timing.yml
```

For a complete output set produced before this runmanager existed, create the
same strict state record while validating it:

```bash
python3 runmanager/runmanager.py validate \
  runmanager/examples/current_geometry_timing.yml --adopt-existing
```

Adoption is refused unless every geometry x primary output declared by the
manifest is present.

Validation requires a non-zombie, non-recovered ROOT file and checks:

- exactly the manifest event count in `evt`, `calarr`, `th`, and `tlc`;
- exactly one `runmeta` entry;
- exact branch schema for `evt`, `calhit`, `target`, `calarr`, `runmeta`, `th`,
  and `tlc`;
- BGO vector length `n_layer*n_sector` and TH/TLC vector length 30;
- manifest geometry mode, segmentation, photon-counter mode, primary, seed,
  and output path;
- for v6, beam species/momentum, overlay flag, target material/areal density/
  density/derived length, and π⁻/π⁰ generator momenta;
- for v7, realized-multiplicity branch plus exact multiplicity, gate, and
  profile metadata;
- fixed `physicsFlag=4`, `neutronScale=2`, `inelasticBias=3`, and
  `pionInelasticXSScale=1.65`.

Published BGOegg frusta use `geometry_mode: bgoegg_frustum` and
`segmentation: bgoegg_published`.  Only the published 22-ring layout and the
segment-size-preserving 31-ring extension are accepted.  The latter adds five
forward and four backward rings while keeping 60 azimuthal sectors.

State and LSF logs live under the gitignored
`beta_org/runmanager/.state/`. State includes the manifest SHA-256, job IDs,
LSF state, output path, validation result, beta repository commit, repository
dirty flag, build executable path, and build executable SHA-256. Reusing a
state after the manifest, git commit/dirty state, build path, or binary hash
changes is rejected. Use a new tag for a rebuilt executable.

`runmeta` directly verifies `physicsFlag=4`, `neutronScale=2`,
`inelasticBias=3`, and `pionInelasticXSScale=1.65`. The low-energy neutron
factor 0.5 and nuclear-like gamma factor 0.5 are not separate `runmeta`
columns. They are fixed by the recorded source commit/dirty flag and exact
build executable SHA-256; changing/rebuilding that implementation requires a
new manifest tag and state.

After `status` and/or `validate`, only failed or missing jobs may be submitted
again:

```bash
python3 runmanager/runmanager.py submit \
  runmanager/examples/current_geometry_timing.yml --retry
```

Pending, running, adopted-valid, and valid jobs are never resubmitted. To change a manifest,
use a new tag so output/state from different schemas or configurations cannot
be mixed.
