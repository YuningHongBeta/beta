#!/usr/bin/env bash
set -euo pipefail

task_build=/gpfs/home/had/yhong/geant4/beta/beta_org/build_bgoegg_frustum
task_macro=/gpfs/home/had/yhong/geant4/beta/beta_org/macros/run.mac
task_outdir=${task_build}/output/scan
task_seed=7232026

mkdir -p "${task_outdir}"
cd "${task_build}"

submit_one() {
  local geometry=$1
  local primary=$2
  local tag=$3
  local n_layer=$4
  local n_sector=$5
  local segmentation=$6
  local bgo_offset_cm=$7
  local pc_front_cm=$8
  local pc_outer_deg=$9

  local output_base=${task_outdir}/${tag}
  if [[ -f ${output_base}.root ]]; then
    echo "SKIP ${tag}: ROOT already exists"
    return
  fi

  bsub -J "${tag}" -n 2 -R "span[hosts=1]" -q s -W 120 \
    -oo "${output_base}.out" -eo "${output_base}.err" \
    env \
      BETA_GEOMETRY="${geometry}" \
      BETA_N_LAYER="${n_layer}" \
      BETA_N_SECTOR="${n_sector}" \
      BETA_SEGMENTATION="${segmentation}" \
      BETA_BGO_Z_OFFSET_CM="${bgo_offset_cm}" \
      BETA_PHOTON_COUNTER=upstream \
      BETA_PC_N_LAYERS=3 \
      BETA_PC_PB_THICKNESS_MM=6 \
      BETA_PC_SCINTI_THICKNESS_MM=4 \
      BETA_PC_Z_FRONT_CM="${pc_front_cm}" \
      BETA_PC_DOWN_THETA_INNER_DEG=0.01 \
      BETA_PC_DOWN_THETA_OUTER_DEG=12 \
      BETA_PC_UP_THETA_INNER_DEG=0.01 \
      BETA_PC_UP_THETA_OUTER_DEG="${pc_outer_deg}" \
      BETA_PC_SQUARE_HOLE_MM=60 \
      BETA_PRIMARY="${primary}" \
      BETA_THREADS=2 \
      BETA_SEED="${task_seed}" \
      BETA_WRITE_CALHIT=0 \
      BETA_OUTPUT="output/scan/${tag}" \
      "${task_build}/beta" "${task_macro}"
}

for task_primary in e pim pi0; do
  submit_one current "${task_primary}" \
    "bgoc_balanced_pc3_h60_s7232026_${task_primary}" \
    15 15 uniform_theta 0 55.5 24
  submit_one bgoegg_frustum "${task_primary}" \
    "bgoegg31_balanced_pc3_h60_s7232026_${task_primary}" \
    31 60 bgoegg_published -10 56.0 14
done
