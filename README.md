# beta

Geant4 simulation for a feasibility study of measuring Lambda's beta-decay rate inside nuclei.

The detector concept consists of a spherical BGO calorimeter (R=30-50 cm, 15x15 segmentation) surrounding an enriched ⁶Li target, with plastic scintillator hodoscopes (TH, TLC) for charged-particle tracking and beam-hole openings for upstream/downstream spectrometers.

## Contents

| Directory | Description |
|-----------|-------------|
| [beta_org/](beta_org/) | Current version (Geant4 v11): upgraded physics lists (INCLXX, biasing), selectable primaries (Λ beta decay / π⁻ / π⁰) |
| [Cal_src/](Cal_src/) | Original version (Geant4 v10): initial detector response study with the same geometry |

## Local build and ROOT-output convention

The repository-level `build/` directory is a pre-reorganization legacy build.
Its executable predates the configurable BGO geometry and has no `runmeta` tree;
do not use it for current BGOegg or segmentation studies.  Configure current
builds with `cmake -S beta_org -B beta_org/<build-name>`.

On the KEK `yhong` workspace, bulk ROOT output must be placed on group storage
through

```text
build/root -> /group/had/sks/Users/yhong/rootfile
```

Existing production build `output` directories are symlinks into
`build/root/beta/`.  When creating a new build directory, create its `output`
as a symlink into a uniquely named directory below `build/root/beta/` before
submitting jobs.  Keep only small logs, JSON, and temporary analysis products
in HOME.
