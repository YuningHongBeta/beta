#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 5 && $# -ne 7 && $# -ne 8 ]]; then
  echo "usage: $0 TAG N_LAYER N_SECTOR SEGMENTATION PRIMARY [GEOMETRY PHOTON_COUNTER [BGO_Z_OFFSET_CM]]" >&2
  exit 2
fi

tag=$1
n_layer=$2
n_sector=$3
segmentation=$4
primary=$5
geometry=${6:-current}
photon_counter=${7:-none}
bgo_z_offset_cm=${8:-}

project_dir=$(cd "$(dirname "$0")/.." && pwd)
build_dir=${BETA_BUILD_DIR:-"$project_dir/build_opt"}
if [[ ! -x "$build_dir/beta" ]]; then
  echo "beta executable not found: $build_dir/beta" >&2
  exit 3
fi
mkdir -p "$build_dir/output/scan"

export BETA_N_LAYER="$n_layer"
export BETA_N_SECTOR="$n_sector"
export BETA_SEGMENTATION="$segmentation"
export BETA_PRIMARY="$primary"
export BETA_GEOMETRY="$geometry"
export BETA_PHOTON_COUNTER="$photon_counter"
if [[ -n "$bgo_z_offset_cm" ]]; then
  export BETA_BGO_Z_OFFSET_CM="$bgo_z_offset_cm"
else
  unset BETA_BGO_Z_OFFSET_CM
fi
export BETA_OUTPUT="output/scan/${tag}_${primary}"
export BETA_WRITE_CALHIT=0
export BETA_THREADS="${BETA_THREADS:-8}"
export BETA_SEED="${BETA_SEED:-6302026}"

cd "$build_dir"
exec ./beta macros/run.mac
