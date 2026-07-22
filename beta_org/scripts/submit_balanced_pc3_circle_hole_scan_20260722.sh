#!/usr/bin/env bash
set -euo pipefail

task_build=/gpfs/home/had/yhong/geant4/beta/beta_org/build_bgoegg_frustum
task_macro=/gpfs/home/had/yhong/geant4/beta/beta_org/macros/run.mac
task_outdir=${task_build}/output/scan

mkdir -p "${task_outdir}"
cd "${task_build}"

submit_one() {
  local seed=$1
  local hole_mm=$2
  local geometry=$3
  local primary=$4
  local tag=$5
  local n_layer=$6
  local n_sector=$7
  local segmentation=$8
  local bgo_offset_cm=$9
  local pc_front_cm=${10}
  local pc_outer_deg=${11}

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
      BETA_PC_CIRCULAR_HOLE_MM="${hole_mm}" \
      BETA_PRIMARY="${primary}" \
      BETA_THREADS=2 \
      BETA_SEED="${seed}" \
      BETA_WRITE_CALHIT=0 \
      BETA_OUTPUT="output/scan/${tag}" \
      "${task_build}/beta" "${task_macro}"
}

for task_seed in 7232026 7242026; do
  for task_hole in 60 100 150; do
    for task_primary in e pim pi0; do
      submit_one "${task_seed}" "${task_hole}" current "${task_primary}" \
        "bgoc_balanced_pc3_hc${task_hole}_s${task_seed}_${task_primary}" \
        15 15 uniform_theta 0 55.5 24
      submit_one "${task_seed}" "${task_hole}" bgoegg_frustum "${task_primary}" \
        "bgoegg31_balanced_pc3_hc${task_hole}_s${task_seed}_${task_primary}" \
        31 60 bgoegg_published -10 56.0 14
    done
  done
done
