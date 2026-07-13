#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 5 ]]; then
  echo "usage: $0 TAG N_LAYER N_SECTOR SEGMENTATION PRIMARY" >&2
  exit 2
fi

tag=$1
n_layer=$2
n_sector=$3
segmentation=$4
primary=$5

project_dir=$(cd "$(dirname "$0")/.." && pwd)
build_dir="$project_dir/build_opt"
mkdir -p "$build_dir/output/scan"

export BETA_N_LAYER="$n_layer"
export BETA_N_SECTOR="$n_sector"
export BETA_SEGMENTATION="$segmentation"
export BETA_PRIMARY="$primary"
export BETA_OUTPUT="output/scan/${tag}_${primary}"
export BETA_WRITE_CALHIT=0
export BETA_THREADS="${BETA_THREADS:-8}"
export BETA_SEED="${BETA_SEED:-6302026}"

cd "$build_dir"
exec ./beta macros/run.mac