#!/usr/bin/env python3
"""Freeze and confirm the stand-alone PC point with a pi0-only sample."""

from __future__ import annotations

import argparse
import json
import math
import subprocess
from pathlib import Path

import numpy as np
import yaml

from bgo_features_v2 import load_bgo2


PROJECT = Path(__file__).resolve().parent.parent
THRESHOLDS_MEV = (0.0, 0.1, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0)


def wilson(numerator: int, denominator: int) -> list[float] | None:
    if denominator <= 0:
        return None
    p = numerator / denominator
    z = 1.959963984540054
    scale = 1.0 + z * z / denominator
    center = (p + z * z / (2.0 * denominator)) / scale
    half = z * math.sqrt(
        p * (1.0 - p) / denominator + z * z / (4.0 * denominator**2)
    ) / scale
    return [float(center - half), float(center + half)]


def measured(mask: np.ndarray, condition: np.ndarray) -> dict:
    denominator = int(np.count_nonzero(condition))
    numerator = int(np.count_nonzero(mask & condition))
    return {
        "value": None if denominator == 0 else numerator / denominator,
        "numerator": numerator,
        "denominator": denominator,
        "wilson_95ci": wilson(numerator, denominator),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--scan-json", type=Path, required=True)
    parser.add_argument("--geometry")
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--extractor", type=Path,
                        default=PROJECT / "analysis/bgo_extract_features_v2")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    manifest = yaml.safe_load(args.manifest.read_text(encoding="utf-8"))
    scan = json.loads(args.scan_json.read_text(encoding="utf-8"))
    selected = scan["selected_ideal_upper_bound"]
    design = manifest["geometries"][0]
    frozen_geometry = args.geometry or selected["geometry"]
    if frozen_geometry not in scan["geometries"]:
        raise ValueError("confirmation geometry is absent from screening scan")
    if frozen_geometry != design["name"]:
        raise ValueError("confirmation geometry differs from frozen scan selection")
    if manifest["primaries"] != ["pi0"]:
        raise ValueError("confirmation manifest must be pi0-only")
    if not args.features.is_file():
        args.features.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run([
            str(args.extractor.resolve()), str(args.root.resolve()),
            str(args.features.resolve()), "3",
        ], check=True)

    data, meta, feature_manifest = load_bgo2(args.features, copy=False)
    if meta["nrow"] != manifest["events"]:
        raise ValueError("confirmation row count differs from manifest")
    if not math.isclose(float(meta["threshold_MeV"]), 3.0):
        raise ValueError("BGO threshold is not 3 MeV")
    col = {name: index for index, name in enumerate(feature_manifest["features"])}
    required = {
        "nCl4", "thGeomPath_mm", "thTotalE_MeV",
        "thMatchedDEdx_MeV_per_mm", "tlcAnySegHitGe1_eff05",
        "absDeltaZValid", "absDeltaZLt90", "pcSumE_MeV", "pcGammaN",
        "pcGammaDownN", "pcGammaUpN",
    }
    missing = sorted(required - set(col))
    if missing:
        raise ValueError(f"missing confirmation features: {missing}")

    all_events = np.ones(len(data), dtype=bool)
    bgo = data[:, col["nCl4"]] == 1
    th = (
        (data[:, col["thGeomPath_mm"]] > 0.0)
        & (data[:, col["thTotalE_MeV"]] >= 0.7)
        & (data[:, col["thTotalE_MeV"]] <= 3.0)
        & (data[:, col["thMatchedDEdx_MeV_per_mm"]] >= 0.12)
        & (data[:, col["thMatchedDEdx_MeV_per_mm"]] <= 0.25)
    )
    tlc = data[:, col["tlcAnySegHitGe1_eff05"]] > 0.5
    dz = (
        (data[:, col["absDeltaZValid"]] > 0.5)
        & (data[:, col["absDeltaZLt90"]] > 0.5)
    )
    conditions = {
        "all_generated": all_events,
        "given_bgo_ncluster1": bgo,
        "given_full_pre_pc": bgo & th & tlc & dz,
    }
    gamma = data[:, col["pcGammaN"]] > 0.0
    gamma_down = data[:, col["pcGammaDownN"]] > 0.0
    gamma_up = data[:, col["pcGammaUpN"]] > 0.0
    rows = []
    for threshold in THRESHOLDS_MEV:
        hit = (
            data[:, col["pcSumE_MeV"]] > 0.0 if threshold == 0.0
            else data[:, col["pcSumE_MeV"]] >= threshold
        )
        rows.append({
            "pc_plastic_threshold_MeV": threshold,
            "conditions": {
                label: {
                    "gamma_intercept": measured(gamma, condition),
                    "effective_pc_hit": measured(hit, condition),
                    "response_given_gamma_intercept": measured(
                        hit, condition & gamma),
                    "hit_without_recorded_gamma_intercept": measured(
                        hit & ~gamma, condition),
                }
                for label, condition in conditions.items()
            },
        })

    fixed_threshold = float(scan["screen_threshold_MeV"])
    fixed = next(row for row in rows
                 if row["pc_plastic_threshold_MeV"] == fixed_threshold)
    result = {
        "status": "complete",
        "analysis": "bgoegg_pc_confirm_evaluate",
        "confirmation_policy": (
            "geometry and PC threshold frozen from seed-650 screening; "
            "stand-alone hard PC hit veto; no BGO+PC classifier"
        ),
        "tag": manifest["tag"],
        "seed": manifest["seed"],
        "events": manifest["events"],
        "bgo_threshold_MeV": 3.0,
        "frozen_geometry": design,
        "frozen_pc_plastic_threshold_MeV": fixed_threshold,
        "selection_counts": {
            label: int(np.count_nonzero(condition))
            for label, condition in conditions.items()
        },
        "fixed_point": fixed,
        "operating_points": rows,
        "gamma_intercept_sides_all_generated": {
            "any": measured(gamma, all_events),
            "upstream": measured(gamma_up, all_events),
            "downstream": measured(gamma_down, all_events),
        },
        "source_root": str(args.root.resolve()),
        "source_features": str(args.features.resolve()),
        "scan_json": str(args.scan_json.resolve()),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "selection_counts": result["selection_counts"],
        "fixed_point": fixed["conditions"],
    }, indent=2))
    print(args.output)


if __name__ == "__main__":
    main()
