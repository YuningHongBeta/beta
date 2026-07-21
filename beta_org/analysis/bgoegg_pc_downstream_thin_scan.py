#!/usr/bin/env python3
"""Screen thin downstream-only Photon Counter candidates.

The screening denominator is the raw 3 MeV, four-neighbour Ncluster == 1
pi0 sample.  The much smaller full pre-PC denominator is also reported but is
not used to choose a geometry.  A selected point must be confirmed with an
independent high-statistics sample after TH, TLC, and delta-z.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import yaml

from bgo_features_v2 import load_bgo2
from bgoegg_pc_aperture_evaluate import ensure_features, masks, measured


PI0_BEFORE_PC_PERCENT = 5.553333333333336
PIM_PERCENT = 1.0308000000000002
TOTAL_TARGET_PERCENT = 4.0
PC_THRESHOLD_MEV = 0.5
PC_THRESHOLD_SCAN_MEV = (0.0, 0.05, 0.1, 0.2, 0.5, 1.0)
BGOC_PI0_BGO_SURVIVAL = 0.01666
PI0_OTHER_FACTORS_TO_PERCENT = 100.0 * 0.2 * 0.1 * 0.1 * 0.8 / 0.00048
PB_X0_MM = 5.612
PLASTIC_X0_MM = 424.3


def parse_args() -> argparse.Namespace:
    project = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument(
        "--root-dir", type=Path,
        default=project / "build_bgoegg_frustum" / "output" / "scan",
    )
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument(
        "--extractor", type=Path,
        default=project / "analysis" / "bgo_extract_features_v2",
    )
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def evaluate(args: argparse.Namespace, tag: str, design: dict) -> dict:
    name = str(design["name"])
    root = args.root_dir / f"{tag}_{name}_pi0.root"
    feature = args.cache / name / "pi0.bgo2"
    ensure_features(root.resolve(), feature.resolve(), args.extractor.resolve())
    data, metadata, manifest = load_bgo2(feature, copy=False)
    if not math.isclose(float(metadata["threshold_MeV"]), 3.0):
        raise ValueError(f"{name}: BGO threshold is not 3 MeV")
    names = list(manifest["features"])
    columns = {feature_name: index for index, feature_name in enumerate(names)}
    selection = masks(data, names)
    pc_mode = str(design["photon_counter"])
    if pc_mode == "downstream":
        side = "Down"
        outer_key = "pc_down_theta_outer_deg"
    elif pc_mode == "upstream":
        side = "Up"
        outer_key = "pc_up_theta_outer_deg"
    else:
        raise ValueError(f"{name}: expected one-sided PC, got {pc_mode}")
    intercept = data[:, columns[f"pcGamma{side}N"]] > 0.0
    pc_energy = data[:, columns[f"pc{side}E_MeV"]]
    pc_hit = pc_energy >= PC_THRESHOLD_MEV

    bgo_intercept = measured(intercept, selection["bgo"])
    bgo_hit = measured(pc_hit, selection["bgo"])
    bgo_response = measured(pc_hit, selection["bgo"] & intercept)
    pre_intercept = measured(intercept, selection["pre_pc"])
    pre_hit = measured(pc_hit, selection["pre_pc"])
    veto = float(bgo_hit["value"])
    n_generated = int(data.shape[0])
    joint_count = int(np.count_nonzero(selection["bgo"] & ~pc_hit))
    joint_survival = joint_count / n_generated
    joint_measurement = measured(
        selection["bgo"] & ~pc_hit,
        np.ones(n_generated, dtype=bool),
    )
    direct_pi0_percent = PI0_OTHER_FACTORS_TO_PERCENT * joint_survival
    direct_total_percent = direct_pi0_percent + PIM_PERCENT
    direct_effective_veto = 1.0 - joint_survival / BGOC_PI0_BGO_SURVIVAL
    threshold_scan = []
    for threshold in PC_THRESHOLD_SCAN_MEV:
        threshold_hit = (
            pc_energy > 0.0 if threshold == 0.0
            else pc_energy >= threshold
        )
        threshold_joint = measured(
            selection["bgo"] & ~threshold_hit,
            np.ones(n_generated, dtype=bool),
        )
        threshold_scan.append({
            "threshold_MeV": threshold,
            "bgo_and_pc_no_hit_count": threshold_joint["numerator"],
            "pi_total_to_beta_percent": (
                PI0_OTHER_FACTORS_TO_PERCENT * threshold_joint["value"]
                + PIM_PERCENT
            ),
            "pi_total_to_beta_percent_pc_mc_stat_only_95ci": [
                PI0_OTHER_FACTORS_TO_PERCENT * edge + PIM_PERCENT
                for edge in threshold_joint["wilson_95ci"]
            ],
        })

    layers = int(design["pc_n_layers"])
    pb_mm = float(design["pc_pb_thickness_mm"])
    plastic_mm = float(design["pc_scinti_thickness_mm"])
    z_front_mm = 10.0 * float(design["pc_z_front_cm"])
    stack_mm = layers * (pb_mm + plastic_mm)
    z_back_mm = z_front_mm + stack_mm
    outer_deg = float(design[outer_key])
    outer_radius_mm = z_back_mm * math.tan(math.radians(outer_deg))
    material_x0 = layers * (
        pb_mm / PB_X0_MM + plastic_mm / PLASTIC_X0_MM
    )
    total_percent = PI0_BEFORE_PC_PERCENT * (1.0 - veto) + PIM_PERCENT
    return {
        "geometry": name,
        "beta_pc_mode": pc_mode,
        "design": design,
        "source_root": str(root.resolve()),
        "source_features": str(feature.resolve()),
        "bgo_ncluster1_count": int(np.count_nonzero(selection["bgo"])),
        "full_pre_pc_count": int(np.count_nonzero(selection["pre_pc"])),
        "bgo_conditional": {
            "gamma_intercept": bgo_intercept,
            "effective_pc_veto": bgo_hit,
            "response_given_intercept": bgo_response,
        },
        "full_pre_pc_diagnostic": {
            "gamma_intercept": pre_intercept,
            "effective_pc_veto": pre_hit,
        },
        "material": {
            "stack_thickness_mm": stack_mm,
            "total_pb_mm": layers * pb_mm,
            "total_plastic_mm": layers * plastic_mm,
            "approx_radiation_lengths_if_traversed": material_x0,
            "outer_radius_at_back_mm": outer_radius_mm,
        },
        "bgoc_transfer_screening": {
            "pi0_to_beta_percent": PI0_BEFORE_PC_PERCENT * (1.0 - veto),
            "pim_to_beta_percent": PIM_PERCENT,
            "pi_total_to_beta_percent": total_percent,
            "passes_total_below_4_percent": total_percent < TOTAL_TARGET_PERCENT,
        },
        "bgoc_direct_screening": {
            "generated_count": n_generated,
            "bgo_and_pc_no_hit_count": joint_count,
            "joint_survival": joint_survival,
            "joint_survival_measurement": joint_measurement,
            "effective_pc_veto_vs_no_pc_bgoc": direct_effective_veto,
            "pi0_to_beta_percent": direct_pi0_percent,
            "pi0_to_beta_percent_pc_mc_stat_only_95ci": [
                PI0_OTHER_FACTORS_TO_PERCENT * edge
                for edge in joint_measurement["wilson_95ci"]
            ],
            "pim_to_beta_percent": PIM_PERCENT,
            "pi_total_to_beta_percent": direct_total_percent,
            "pi_total_to_beta_percent_pc_mc_stat_only_95ci": [
                PI0_OTHER_FACTORS_TO_PERCENT * edge + PIM_PERCENT
                for edge in joint_measurement["wilson_95ci"]
            ],
            "passes_total_below_4_percent": direct_total_percent < TOTAL_TARGET_PERCENT,
            "threshold_scan": threshold_scan,
        },
    }


def main() -> None:
    args = parse_args()
    manifest = yaml.safe_load(args.manifest.read_text(encoding="utf-8"))
    tag = str(manifest["tag"])
    geometry_modes = sorted({
        str(item.get("geometry_mode", "current"))
        for item in manifest["geometries"]
    })
    if len(geometry_modes) != 1:
        raise ValueError(f"mixed geometry modes are unsupported: {geometry_modes}")
    geometry_mode = geometry_modes[0]
    rows = [evaluate(args, tag, item) for item in manifest["geometries"]]
    performance_key = (
        "bgoc_direct_screening" if geometry_mode == "current"
        else "bgoc_transfer_screening"
    )
    installable = [
        row for row in rows
        if "reference" not in row["geometry"]
        and row[performance_key]["passes_total_below_4_percent"]
    ]
    selected = min(
        installable,
        key=lambda row: (
            row["material"]["stack_thickness_mm"],
            row["material"]["approx_radiation_lengths_if_traversed"],
            row[performance_key]["pi_total_to_beta_percent"],
        ),
    ) if installable else None
    required_veto = 1.0 - (
        TOTAL_TARGET_PERCENT - PIM_PERCENT
    ) / PI0_BEFORE_PC_PERCENT
    output = {
        "analysis": "thin downstream-only Photon Counter screening",
        "selection": (
            f"geometry_mode={geometry_mode}; 3 MeV/cell; four-neighbour "
            "Ncluster == 1; stand-alone downstream plastic Edep >= 0.5 MeV veto"
        ),
        "target": {
            "pi_total_to_beta_percent": TOTAL_TARGET_PERCENT,
            "required_pc_veto_given_current_bgoc": required_veto,
        },
        "selection_policy": (
            "choose the thinnest non-reference point below 4% in BGO-conditional "
            "screening; for current BGOC use the directly correlated BGO-and-PC "
            "survival; require independent high-statistics confirmation"
        ),
        "selected_screening_candidate": selected,
        "candidates": rows,
        "scope": {
            "screening_only": True,
            "not_included": (
                "K18 accepted-track material crossing, momentum-resolution loss, "
                "supports/readout, rate, pile-up, and BGOC end-to-end correlation"
            ),
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    summary = [{
        "geometry": row["geometry"],
        "thickness_mm": row["material"]["stack_thickness_mm"],
        "x0": row["material"]["approx_radiation_lengths_if_traversed"],
        "intercept": row["bgo_conditional"]["gamma_intercept"]["value"],
        "veto": row["bgo_conditional"]["effective_pc_veto"]["value"],
        "total_percent": row[performance_key]["pi_total_to_beta_percent"],
    } for row in rows]
    print(json.dumps(summary, indent=2))
    print(json.dumps({"selected": None if selected is None else selected["geometry"]}))
    print(args.output)


if __name__ == "__main__":
    main()
