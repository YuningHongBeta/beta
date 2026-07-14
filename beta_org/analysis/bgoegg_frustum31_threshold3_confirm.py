#!/usr/bin/env python3
"""Apply the single selected exact-BGOegg classifier to an independent seed."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from bgo_compare_v1 import (
    FEATURE_GROUPS,
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
    GUARDED_TARGET,
    HURDLE_FEATURE_SETS,
    SPECIES_FILE,
    TARGET,
    fit_hurdle,
    hurdle_feature_values,
    hurdle_scores,
)
from bgoegg_frustum31_threshold3_scan import ROBUST_TARGET


def paths(directory: Path) -> dict[str, Path]:
    result = {key: directory / value for key, value in SPECIES_FILE.items()}
    missing = [str(path) for path in result.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing BGO2 input: " + ", ".join(missing))
    return result


def standard_train_apply(
    training_directory: Path,
    confirmation_directory: Path,
    features: list[str],
    method: str,
) -> dict:
    training_paths = paths(training_directory)
    training = load_geometry(
        training_paths["electron"], training_paths["pim"], training_paths["pi0"],
        "training", features,
    )
    fit_parts = {}
    validation_scores = {}
    for species, sample in training["samples"].items():
        masks = split_masks(sample["event_id"])
        fit_parts[species] = sample["x"][masks["fit"]]
    model = fit_method("qda" if method == "pairwise_qda" else method, fit_parts)
    for species, sample in training["samples"].items():
        mask = split_masks(sample["event_id"])["validation"]
        validation_scores[species] = score_method(method, model, sample["x"][mask])

    confirmation_paths = paths(confirmation_directory)
    confirmation = load_geometry(
        confirmation_paths["electron"], confirmation_paths["pim"],
        confirmation_paths["pi0"], "confirmation", features,
    )
    confirmation_scores = {
        species: score_method(method, model, sample["x"])
        for species, sample in confirmation["samples"].items()
    }
    if method == "pairwise_qda":
        cuts = choose_pairwise_thresholds(validation_scores, *TARGET)
        validation_result = metrics_pairwise(validation_scores, cuts)
        confirmation_result = metrics_pairwise(confirmation_scores, cuts)
    else:
        cuts = choose_threshold(validation_scores, *TARGET)
        validation_result = metrics(validation_scores, cuts)
        confirmation_result = metrics(confirmation_scores, cuts)
    return {
        "classifier_inputs": features,
        "method": method,
        "target": {"pim_reject": TARGET[0], "pi0_reject": TARGET[1]},
        "threshold": (
            float(cuts) if np.isscalar(cuts)
            else [float(value) for value in cuts]
        ),
        "training_validation": validation_result,
        "confirmation_all_100k": confirmation_result,
    }


def hurdle_train_apply(
    training_directory: Path,
    confirmation_directory: Path,
    feature_set: str,
    target: tuple[float, float],
) -> dict:
    features = HURDLE_FEATURE_SETS[feature_set]
    training = hurdle_feature_values(training_directory, "training", features)
    model = fit_hurdle(training)
    validation_scores = {
        species: hurdle_scores(
            model, sample, split_masks(sample["event_id"])["validation"]
        )
        for species, sample in training["samples"].items()
    }
    cuts = choose_pairwise_thresholds(validation_scores, *target)

    confirmation = hurdle_feature_values(
        confirmation_directory, "confirmation", features
    )
    confirmation_scores = {
        species: hurdle_scores(
            model, sample, np.ones(len(sample["event_id"]), dtype=bool)
        )
        for species, sample in confirmation["samples"].items()
    }
    return {
        "classifier_inputs": features,
        "method": "nHit-binned diagonal Gaussian likelihood",
        "target": {"pim_reject": target[0], "pi0_reject": target[1]},
        "threshold": [float(value) for value in cuts],
        "training_validation": metrics_pairwise(validation_scores, cuts),
        "confirmation_all_100k": metrics_pairwise(confirmation_scores, cuts),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scan-json", type=Path, required=True)
    parser.add_argument("--training", type=Path, required=True)
    parser.add_argument("--confirmation", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    scan = json.loads(args.scan_json.read_text(encoding="utf-8"))
    selected = scan["selected"]
    if selected["family"] in {"hurdle", "hurdle_guarded", "hurdle_robust"}:
        target = {
            "hurdle": TARGET,
            "hurdle_guarded": GUARDED_TARGET,
            "hurdle_robust": ROBUST_TARGET,
        }[selected["family"]]
        result = hurdle_train_apply(
            args.training,
            args.confirmation,
            selected["feature_set"],
            target,
        )
    else:
        features = (
            FEATURE_GROUPS["bgo"]
            if selected["family"] == "standard_full"
            else COMPACT_FEATURE_SETS[selected["feature_set"]]
        )
        result = standard_train_apply(
            args.training, args.confirmation, features, selected["method"]
        )
    output = {
        "status": "complete",
        "analysis": "bgoegg_frustum31_threshold3_confirm",
        "confirmation_policy": (
            "geometry, feature set, model family, fit rule, and validation cut are "
            "frozen from seed 6302026; seed 7302026 is evaluated once without refitting"
        ),
        "selected": selected,
        "training_directory": str(args.training.resolve()),
        "confirmation_directory": str(args.confirmation.resolve()),
        **result,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output["confirmation_all_100k"], indent=2))
    print(args.output)


if __name__ == "__main__":
    main()
