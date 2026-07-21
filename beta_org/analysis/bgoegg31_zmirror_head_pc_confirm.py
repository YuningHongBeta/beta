#!/usr/bin/env python3
"""High-statistics confirmation of the selected z-mirrored BGOegg31 PC."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import yaml

from bgo_features_v2 import load_bgo2
from bgoegg_pc_aperture_evaluate import ensure_features, masks, measured


PI0_FACTOR_TO_PERCENT = 100.0 * 0.2 * 0.1 * 0.1 * 0.8 / 0.00048
THRESHOLDS_MEV = (0.1, 0.2, 0.5)


def parse_args() -> argparse.Namespace:
    project = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser()
    parser.add_argument("--screening-json", type=Path, required=True)
    parser.add_argument("--highstat-manifest", type=Path, required=True)
    parser.add_argument("--electron-manifest", type=Path, required=True)
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


def feature_for(manifest_path: Path, species: str, args: argparse.Namespace):
    manifest = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
    design = manifest["geometries"][0]
    tag = str(manifest["tag"])
    name = str(design["name"])
    root = args.root_dir / f"{tag}_{name}_{species}.root"
    feature = args.cache / tag / f"{species}.bgo2"
    ensure_features(root.resolve(), feature.resolve(), args.extractor.resolve())
    data, metadata, feature_manifest = load_bgo2(feature, copy=False)
    if not math.isclose(float(metadata["threshold_MeV"]), 3.0):
        raise ValueError(f"{feature}: expected 3 MeV BGO threshold")
    return data, list(feature_manifest["features"]), design, root


def main() -> None:
    args = parse_args()
    screening = json.loads(args.screening_json.read_text(encoding="utf-8"))
    pim_percent = float(screening["baseline"]["pim_to_beta_percent"])

    pi0, names, design, pi0_root = feature_for(
        args.highstat_manifest, "pi0", args
    )
    columns = {name: index for index, name in enumerate(names)}
    pi0_bgo = masks(pi0, names)["bgo"]
    pi0_pc = pi0[:, columns["pcUpE_MeV"]]
    threshold_scan = []
    for threshold in THRESHOLDS_MEV:
        joint = measured(pi0_bgo & (pi0_pc < threshold))
        total = PI0_FACTOR_TO_PERCENT * joint["value"] + pim_percent
        interval = [
            PI0_FACTOR_TO_PERCENT * edge + pim_percent
            for edge in joint["wilson_95ci"]
        ]
        threshold_scan.append({
            "threshold_MeV": threshold,
            "bgo_and_pc_no_hit": joint,
            "pi0_to_beta_percent": PI0_FACTOR_TO_PERCENT * joint["value"],
            "pim_to_beta_percent": pim_percent,
            "pi_total_to_beta_percent": total,
            "pi_total_to_beta_percent_pc_mc_stat_only_95ci": interval,
        })

    electron, electron_names, _, electron_root = feature_for(
        args.electron_manifest, "e", args
    )
    electron_columns = {
        name: index for index, name in enumerate(electron_names)
    }
    electron_bgo = masks(electron, electron_names)["bgo"]
    electron_pc = electron[:, electron_columns["pcUpE_MeV"]]
    electron_rows = []
    for threshold in THRESHOLDS_MEV:
        joint = measured(electron_bgo & (electron_pc < threshold))
        electron_rows.append({
            "threshold_MeV": threshold,
            "bgo_and_pc_no_hit": joint,
            "pc_false_veto_given_bgo": measured(
                electron_pc >= threshold, electron_bgo
            ),
        })

    layers = int(design["pc_n_layers"])
    pb_mm = float(design["pc_pb_thickness_mm"])
    plastic_mm = float(design["pc_scinti_thickness_mm"])
    z_front_mm = 10.0 * float(design["pc_z_front_cm"])
    stack_mm = layers * (pb_mm + plastic_mm)
    outer_deg = float(design["pc_up_theta_outer_deg"])
    output = {
        "analysis": "BGOegg31 z-mirrored downstream-head PC confirmation",
        "orientation": {
            "k18_upstream_bottom_opening_deg": 5.336032242257286,
            "k18_downstream_head_opening_deg": 12.0,
            "beta_mapping": (
                "K18 z mirror; beta legacy PCUp at -z is physical downstream"
            ),
        },
        "selection": (
            "T3 3 MeV/cell; four-neighbour Ncluster == 1; downstream-head "
            "PC plastic deposited-energy threshold"
        ),
        "design": design,
        "material": {
            "stack_thickness_mm": stack_mm,
            "total_pb_mm": layers * pb_mm,
            "total_plastic_mm": layers * plastic_mm,
            "approx_radiation_lengths_if_traversed": layers * (
                pb_mm / 5.612 + plastic_mm / 424.3
            ),
            "outer_radius_at_back_mm": (
                (z_front_mm + stack_mm) * math.tan(math.radians(outer_deg))
            ),
        },
        "pi0_source_root": str(pi0_root.resolve()),
        "pi0_bgo_ncluster1": measured(pi0_bgo),
        "threshold_scan": threshold_scan,
        "electron_source_root": str(electron_root.resolve()),
        "electron_bgo_ncluster1": measured(electron_bgo),
        "electron_threshold_scan": electron_rows,
        "scope": {
            "pim_term_source": str(args.screening_json.resolve()),
            "reference_factors": (
                "Hong/Kamada TH, TLC, delta-z and branching ratios"
            ),
            "not_included": (
                "optical/electronic response, pile-up, support/readout "
                "material, and K18 momentum-resolution change"
            ),
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "design": design["name"],
        "threshold_scan": threshold_scan,
        "electron_threshold_scan": electron_rows,
    }, indent=2))
    print(args.output)


if __name__ == "__main__":
    main()
