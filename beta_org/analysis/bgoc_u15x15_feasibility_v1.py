#!/usr/bin/env python3
"""Original 15x15 BGOC threshold/cutout feasibility study.

The 1 MeV reference is the unmasked four-neighbour ``nCl4 == 1`` selection.
Three-MeV classifier families are selected with seed-6302026 validation events
and applied once, without refitting, to the independent seed-7302026 sample.
"""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

import numpy as np

from bgo_compare_v1 import (
    FEATURE_GROUPS,
    METHODS,
    analyze_geometry,
    choose_pairwise_thresholds,
    choose_threshold,
    fit_method,
    load_geometry,
    metrics,
    metrics_pairwise,
    score_method,
    split_masks,
)
from bgo_threshold3_analysis_v1 import (
    COMPACT_FEATURE_SETS,
    HURDLE_FEATURE_SETS,
    analyze_hurdle,
    fit_hurdle,
    hurdle_feature_values,
    hurdle_scores,
)


PROJECT = Path(__file__).resolve().parent.parent
SPECIES_FILE = {"electron": "e.bgo2", "pim": "pim.bgo2", "pi0": "pi0.bgo2"}
ROOT_INPUTS = {
    "training": {
        "electron": "incl165_v2_u15x15_e.root",
        "pim": "incl165_v2_u15x15_pim.root",
        "pi0": "bgoc_u15x15_training_s6302026_u15x15_pi0.root",
    },
    "confirmation": {
        "electron": "bgoc_u15x15_t3_confirm_s7302026_u15x15_e.root",
        "pim": "bgoc_u15x15_t3_confirm_s7302026_u15x15_pim.root",
        "pi0": "bgoc_u15x15_t3_confirm_s7302026_u15x15_pi0.root",
    },
}


def triplet(directory: Path) -> dict[str, Path]:
    result = {name: directory / filename for name, filename in SPECIES_FILE.items()}
    missing = [str(path) for path in result.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing BGO2 file: " + ", ".join(missing))
    return result


def geometry(directory: Path, label: str, features: list[str] | None) -> dict:
    paths = triplet(directory)
    return load_geometry(
        paths["electron"], paths["pim"], paths["pi0"], label, features
    )


def mask_sets(overlap_path: Path) -> dict[str, list[int]]:
    overlap = json.loads(overlap_path.read_text(encoding="utf-8"))
    current = overlap["geometries"]["current_bgoc_z0"]
    detectors = current["detectors"]
    bac = set(detectors["BAC_material_only"]["copy_numbers"])
    sdc2 = set(detectors["SDC2_full_gas_box"]["copy_numbers"])
    return {
        "none": [],
        "bac7": sorted(bac),
        "sdc2_42": sorted(sdc2),
        "sdc2_bac49": sorted(sdc2 | bac),
        "sdc1active_bac66": current["union_active_gas_lower_bound"]["copy_numbers"],
        "sdc1full_bac80": current["union_full_gas_box"]["copy_numbers"],
    }


def extract_inputs(args: argparse.Namespace, masks: dict[str, list[int]]) -> None:
    args.cache.mkdir(parents=True, exist_ok=True)
    mask_directory = args.cache / "masks"
    mask_directory.mkdir(exist_ok=True)
    for mask_name, cells in masks.items():
        mask_path = mask_directory / f"{mask_name}.txt"
        mask_path.write_text("".join(f"{cell}\n" for cell in cells), encoding="utf-8")
        selected_samples = args.extract_sample or list(ROOT_INPUTS)
        for sample_name in selected_samples:
            species_roots = ROOT_INPUTS[sample_name]
            for threshold in (1, 3):
                destination = args.cache / "bgo2" / mask_name / f"{sample_name}_T{threshold}"
                destination.mkdir(parents=True, exist_ok=True)
                for species, root_name in species_roots.items():
                    command = [
                        str(args.extractor),
                        str(args.root_directory / root_name),
                        str(destination / SPECIES_FILE[species]),
                        str(threshold),
                        "--allow-legacy-hodo",
                    ]
                    if cells:
                        command.extend(["--mask-cells", str(mask_path)])
                    subprocess.run(command, check=True)


def cluster_metrics(directory: Path, label: str, split: str | None) -> dict:
    loaded = geometry(directory, label, ["nCl4"])
    accepted = {}
    counts = {}
    for species, sample in loaded["samples"].items():
        mask = np.ones(len(sample["event_id"]), dtype=bool)
        if split is not None:
            mask = split_masks(sample["event_id"])[split]
        values = np.asarray(sample["raw"][mask, sample["names"].index("nCl4")])
        accepted[species] = values == 1
        counts[species] = int(len(values))
    return {
        "e_keep": float(accepted["electron"].mean()),
        "pim_reject": float(1.0 - accepted["pim"].mean()),
        "pi0_reject": float(1.0 - accepted["pi0"].mean()),
        "counts": counts,
    }


def selected_standard(
    directory: Path, label: str, features: list[str], target: tuple[float, float]
) -> dict:
    result = analyze_geometry(geometry(directory, label, features), [target], METHODS)
    target_result = next(iter(result["targets"].values()))
    method = target_result["selected_method"]
    return {
        "method": method,
        "validation": target_result["methods"][method]["validation"],
    }


def scan_candidates(
    directory: Path, label: str, target: tuple[float, float]
) -> tuple[dict, list[dict]]:
    rows = []
    full = selected_standard(directory, label + ":full", FEATURE_GROUPS["bgo"], target)
    rows.append({
        "family": "standard_full", "feature_set": "bgo",
        "features": FEATURE_GROUPS["bgo"], **full,
    })
    for feature_set, features in COMPACT_FEATURE_SETS.items():
        result = selected_standard(
            directory, f"{label}:compact:{feature_set}", features, target
        )
        rows.append({
            "family": "compact", "feature_set": feature_set,
            "features": features, **result,
        })
    for feature_set, features in HURDLE_FEATURE_SETS.items():
        result = analyze_hurdle(
            directory, f"{label}:hurdle:{feature_set}", features, target
        )
        rows.append({
            "family": "hurdle", "feature_set": feature_set,
            "features": features,
            "method": "nHit-binned diagonal Gaussian likelihood",
            "validation": result["validation"],
        })
    selected = max(rows, key=lambda row: row["validation"]["e_keep"])
    return selected, rows


def standard_train_apply(
    training_directory: Path,
    confirmation_directory: Path,
    features: list[str],
    method: str,
    target: tuple[float, float],
) -> dict:
    training = geometry(training_directory, "training-selected", features)
    fit_parts = {
        species: sample["x"][split_masks(sample["event_id"])["fit"]]
        for species, sample in training["samples"].items()
    }
    model = fit_method("qda" if method == "pairwise_qda" else method, fit_parts)
    validation_scores = {
        species: score_method(
            method, model, sample["x"][split_masks(sample["event_id"])["validation"]]
        )
        for species, sample in training["samples"].items()
    }
    confirmation = geometry(confirmation_directory, "confirmation-selected", features)
    confirmation_scores = {
        species: score_method(method, model, sample["x"])
        for species, sample in confirmation["samples"].items()
    }
    if method == "pairwise_qda":
        cuts = choose_pairwise_thresholds(validation_scores, *target)
        validation_result = metrics_pairwise(validation_scores, cuts)
        confirmation_result = metrics_pairwise(confirmation_scores, cuts)
    else:
        cuts = choose_threshold(validation_scores, *target)
        validation_result = metrics(validation_scores, cuts)
        confirmation_result = metrics(confirmation_scores, cuts)
    return {
        "target": {"pim_reject": target[0], "pi0_reject": target[1]},
        "threshold": (
            float(cuts) if np.isscalar(cuts) else [float(value) for value in cuts]
        ),
        "training_validation": validation_result,
        "confirmation_all": confirmation_result,
    }


def hurdle_train_apply(
    training_directory: Path,
    confirmation_directory: Path,
    features: list[str],
    target: tuple[float, float],
) -> dict:
    training = hurdle_feature_values(training_directory, "training-hurdle", features)
    model = fit_hurdle(training)
    validation_scores = {
        species: hurdle_scores(
            model, sample, split_masks(sample["event_id"])["validation"]
        )
        for species, sample in training["samples"].items()
    }
    cuts = choose_pairwise_thresholds(validation_scores, *target)
    confirmation = hurdle_feature_values(
        confirmation_directory, "confirmation-hurdle", features
    )
    confirmation_scores = {
        species: hurdle_scores(
            model, sample, np.ones(len(sample["event_id"]), dtype=bool)
        )
        for species, sample in confirmation["samples"].items()
    }
    return {
        "target": {"pim_reject": target[0], "pi0_reject": target[1]},
        "threshold": [float(value) for value in cuts],
        "training_validation": metrics_pairwise(validation_scores, cuts),
        "confirmation_all": metrics_pairwise(confirmation_scores, cuts),
    }


def apply_selected(
    selected: dict,
    training_directory: Path,
    confirmation_directory: Path,
    target: tuple[float, float],
) -> dict:
    if selected["family"] == "hurdle":
        return hurdle_train_apply(
            training_directory, confirmation_directory, selected["features"], target
        )
    return standard_train_apply(
        training_directory,
        confirmation_directory,
        selected["features"],
        selected["method"],
        target,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--cache", type=Path, default=PROJECT / "tmp" / "bgoc_u15x15_feasibility_v1"
    )
    parser.add_argument(
        "--overlap", type=Path,
        default=Path("/gpfs/home/had/yhong/geant4/k18geant4/param/BGOEGG/")
        / "bgo_material_overlap_comparison_e63_20260715.json",
    )
    parser.add_argument(
        "--root-directory", type=Path, default=PROJECT / "build_opt/output/scan"
    )
    parser.add_argument(
        "--extractor", type=Path, default=PROJECT / "analysis/bgo_extract_features_v2.masktest"
    )
    parser.add_argument("--extract", action="store_true")
    parser.add_argument("--extract-only", action="store_true")
    parser.add_argument(
        "--extract-sample", action="append", choices=tuple(ROOT_INPUTS),
        help="limit --extract to training or confirmation; repeatable",
    )
    parser.add_argument(
        "--output", type=Path,
        default=PROJECT / "analysis/results/bgoegg_v1/bgoc_u15x15_feasibility_v1.json",
    )
    args = parser.parse_args()
    args.cache = args.cache.resolve()
    masks = mask_sets(args.overlap)
    if args.extract:
        extract_inputs(args, masks)
        if args.extract_only:
            return

    bgo2 = args.cache / "bgo2"
    reference = cluster_metrics(
        bgo2 / "none/training_T1", "reference-training-T1", "validation"
    )
    target = (reference["pim_reject"], reference["pi0_reject"])
    guarded_target = (min(0.9999, target[0] + 0.005), min(0.9999, target[1] + 0.005))
    results = {}
    for mask_name, cells in masks.items():
        training_t3 = bgo2 / mask_name / "training_T3"
        confirmation_t3 = bgo2 / mask_name / "confirmation_T3"
        selected, candidates = scan_candidates(training_t3, mask_name, target)
        results[mask_name] = {
            "disabled_cells": cells,
            "disabled_cell_count": len(cells),
            "coverage_remaining_fraction": 1.0 - len(cells) / 225.0,
            "ncl4_eq_1": {
                "training_validation_T1": cluster_metrics(
                    bgo2 / mask_name / "training_T1", mask_name + ":train:T1", "validation"
                ),
                "confirmation_T1": cluster_metrics(
                    bgo2 / mask_name / "confirmation_T1", mask_name + ":confirm:T1", None
                ),
                "confirmation_T3": cluster_metrics(
                    confirmation_t3, mask_name + ":confirm:T3", None
                ),
            },
            "selection_rule": (
                "maximum seed-6302026 validation electron retention among declared "
                "standard, compact, and hurdle candidates at the unmasked T1 target"
            ),
            "selected": {
                key: selected[key]
                for key in ("family", "feature_set", "features", "method", "validation")
            },
            "exact_target": apply_selected(
                selected, training_t3, confirmation_t3, target
            ),
            "guarded_target": apply_selected(
                selected, training_t3, confirmation_t3, guarded_target
            ),
            "candidate_table": candidates,
        }
    output = {
        "status": "complete",
        "analysis": "bgoc_u15x15_feasibility_v1",
        "geometry": "original uniform-theta 15x15 BGOC; 225 cells",
        "threshold_policy": (
            "hard per-cell threshold in feature extraction; sub-threshold deposits are absent"
        ),
        "reference_definition": (
            "unmasked 1 MeV, four-neighbour reconstructed cluster count nCl4 == 1; "
            "target rates frozen from seed-6302026 validation events"
        ),
        "reference_training_validation": reference,
        "target": {"pim_reject": target[0], "pi0_reject": target[1]},
        "guarded_target": {
            "pim_reject": guarded_target[0], "pi0_reject": guarded_target[1],
            "definition": "fixed +0.5 percentage-point validation margin before confirmation",
        },
        "split": {
            "fit": "seed 6302026 eventID mod 4 == 0",
            "validation": "seed 6302026 eventID mod 4 == 2; candidate and cut selection",
            "confirmation": "all 100k events from independent seed 7302026; no refit",
        },
        "overlap_source": str(args.overlap.resolve()),
        "results": results,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    concise = {
        name: {
            "disabled": value["disabled_cell_count"],
            "selected": value["selected"]["family"] + ":" + value["selected"]["feature_set"],
            "exact_confirmation": value["exact_target"]["confirmation_all"],
            "guarded_confirmation": value["guarded_target"]["confirmation_all"],
        }
        for name, value in results.items()
    }
    print(json.dumps({"target": output["target"], "results": concise}, indent=2))
    print(args.output)


if __name__ == "__main__":
    main()
