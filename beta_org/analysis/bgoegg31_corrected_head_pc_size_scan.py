#!/usr/bin/env python3
"""Evaluate corrected BGOegg31 head-side Photon Counter outer size."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import yaml

from bgo_features_v2 import load_bgo2
from bgoegg_pc_aperture_evaluate import ensure_features, masks, measured


THRESHOLDS_MEV = (0.1, 0.2, 0.5)
PI0_FACTOR_TO_PERCENT = 100.0 * 0.2 * 0.1 * 0.1 * 0.8 / 0.00048
PIM_FACTOR_TO_PERCENT = 100.0 * 0.4 * 0.001 * 0.045 / 0.00048
PB_X0_MM = 5.612
PLASTIC_X0_MM = 424.3


def parse_args() -> argparse.Namespace:
    project = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline-manifest", type=Path, required=True)
    parser.add_argument("--scan-manifest", type=Path, required=True)
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


def load_feature(root: Path, feature: Path, extractor: Path):
    ensure_features(root.resolve(), feature.resolve(), extractor.resolve())
    data, metadata, manifest = load_bgo2(feature, copy=False)
    if not math.isclose(float(metadata["threshold_MeV"]), 3.0):
        raise ValueError(f"{feature}: expected 3 MeV BGO threshold")
    names = list(manifest["features"])
    return data, names


def main() -> None:
    args = parse_args()
    baseline_manifest = yaml.safe_load(
        args.baseline_manifest.read_text(encoding="utf-8")
    )
    scan_manifest = yaml.safe_load(args.scan_manifest.read_text(encoding="utf-8"))
    baseline_tag = str(baseline_manifest["tag"])
    baseline_name = str(baseline_manifest["geometries"][0]["name"])

    baseline = {}
    for species in ("e", "pim", "pi0"):
        root = args.root_dir / f"{baseline_tag}_{baseline_name}_{species}.root"
        feature = args.cache / "baseline" / f"{species}.bgo2"
        data, names = load_feature(root, feature, args.extractor)
        selection = masks(data, names)["bgo"]
        baseline[species] = measured(selection)

    pim_percent = PIM_FACTOR_TO_PERCENT * baseline["pim"]["value"]
    pi0_no_pc_percent = PI0_FACTOR_TO_PERCENT * baseline["pi0"]["value"]
    scan_tag = str(scan_manifest["tag"])
    candidates = []
    for design in scan_manifest["geometries"]:
        name = str(design["name"])
        root = args.root_dir / f"{scan_tag}_{name}_pi0.root"
        feature = args.cache / name / "pi0.bgo2"
        data, names = load_feature(root, feature, args.extractor)
        columns = {feature_name: index for index, feature_name in enumerate(names)}
        bgo = masks(data, names)["bgo"]
        pc_energy = data[:, columns["pcUpE_MeV"]]
        intercept = data[:, columns["pcGammaUpN"]] > 0.0
        threshold_rows = []
        for threshold in THRESHOLDS_MEV:
            pc_hit = pc_energy >= threshold
            joint = measured(bgo & ~pc_hit)
            total = PI0_FACTOR_TO_PERCENT * joint["value"] + pim_percent
            interval = [
                PI0_FACTOR_TO_PERCENT * edge + pim_percent
                for edge in joint["wilson_95ci"]
            ]
            threshold_rows.append({
                "threshold_MeV": threshold,
                "bgo_and_pc_no_hit": joint,
                "pi0_to_beta_percent": PI0_FACTOR_TO_PERCENT * joint["value"],
                "pim_to_beta_percent": pim_percent,
                "pi_total_to_beta_percent": total,
                "pi_total_to_beta_percent_pc_mc_stat_only_95ci": interval,
            })

        layers = int(design["pc_n_layers"])
        pb_mm = float(design["pc_pb_thickness_mm"])
        plastic_mm = float(design["pc_scinti_thickness_mm"])
        z_front_mm = 10.0 * float(design["pc_z_front_cm"])
        stack_mm = layers * (pb_mm + plastic_mm)
        z_back_mm = z_front_mm + stack_mm
        outer_deg = float(design["pc_up_theta_outer_deg"])
        candidates.append({
            "geometry": name,
            "design": design,
            "source_root": str(root.resolve()),
            "bgo_ncluster1": measured(bgo),
            "gamma_intercept_given_bgo": measured(intercept, bgo),
            "pc_veto_given_bgo_at_0p2_MeV": measured(pc_energy >= 0.2, bgo),
            "threshold_scan": threshold_rows,
            "material": {
                "stack_thickness_mm": stack_mm,
                "total_pb_mm": layers * pb_mm,
                "total_plastic_mm": layers * plastic_mm,
                "approx_radiation_lengths_if_traversed": layers * (
                    pb_mm / PB_X0_MM + plastic_mm / PLASTIC_X0_MM
                ),
                "outer_radius_at_back_mm": (
                    z_back_mm * math.tan(math.radians(outer_deg))
                ),
            },
        })

    def at_0p2(row: dict) -> dict:
        return next(
            point for point in row["threshold_scan"]
            if point["threshold_MeV"] == 0.2
        )

    passing = [
        row for row in candidates
        if at_0p2(row)["pi_total_to_beta_percent"] < 4.0
    ]
    minimum_passing = min(
        passing,
        key=lambda row: (
            float(row["design"]["pc_up_theta_outer_deg"]),
            at_0p2(row)["pi_total_to_beta_percent"],
        ),
    ) if passing else None
    best_total = min(
        at_0p2(row)["pi_total_to_beta_percent"] for row in candidates
    )
    plateau = [
        row for row in passing
        if at_0p2(row)["pi_total_to_beta_percent"] <= best_total + 0.02
    ]
    selected = min(
        plateau,
        key=lambda row: float(row["design"]["pc_up_theta_outer_deg"]),
    ) if plateau else min(
        candidates, key=lambda row: at_0p2(row)["pi_total_to_beta_percent"]
    )
    output = {
        "analysis": "corrected BGOegg31 downstream-head PC outer-size scan",
        "orientation": {
            "beta_upstream_bottom_opening_deg": 5.336032242257286,
            "beta_downstream_head_opening_deg": 12.0,
            "k18_mapping": "z mirror: bottom upstream -z, head downstream +z",
        },
        "selection": (
            "T3 3 MeV/cell; four-neighbour Ncluster == 1; physical downstream "
            "head-side PC; plastic Edep threshold scan"
        ),
        "baseline": {
            "measurements": baseline,
            "pi0_no_pc_to_beta_percent": pi0_no_pc_percent,
            "pim_to_beta_percent": pim_percent,
            "pi_total_no_pc_to_beta_percent": pi0_no_pc_percent + pim_percent,
        },
        "selection_policy": (
            "screening: among central totals below 4% at 0.2 MeV, choose the "
            "smallest PC outer angle within 0.02 percentage point of the best "
            "screening total (the saturated plateau); independent "
            "high-statistics confirmation required"
        ),
        "minimum_below_4_percent_candidate": minimum_passing,
        "selected_screening_candidate": selected,
        "candidates": candidates,
        "scope": {
            "screening_only": True,
            "reference_factors": "Hong/Kamada TH, TLC, delta-z and branching ratios",
            "not_included": (
                "optical/electronic response, pile-up, support/readout material, "
                "and K18 momentum-resolution change"
            ),
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "baseline": output["baseline"],
        "candidates": [{
            "geometry": row["geometry"],
            "outer_deg": row["design"]["pc_up_theta_outer_deg"],
            "outer_radius_mm": row["material"]["outer_radius_at_back_mm"],
            "total_0p2_percent": at_0p2(row)["pi_total_to_beta_percent"],
        } for row in candidates],
        "selected": selected["geometry"],
    }, indent=2))
    print(args.output)


if __name__ == "__main__":
    main()
