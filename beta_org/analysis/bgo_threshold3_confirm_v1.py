#!/usr/bin/env python3
"""Apply frozen 3 MeV hurdle classifiers to independent-seed samples."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from bgo_compare_v1 import choose_pairwise_thresholds, metrics_pairwise, split_masks
from bgo_threshold3_analysis_v1 import (
    GUARDED_TARGET,
    HURDLE_FEATURE_SETS,
    PROJECT,
    TARGET,
    fit_hurdle,
    hurdle_feature_values,
    hurdle_scores,
)


DEFAULT_CACHE = PROJECT / "tmp" / "bgo_threshold3_v1"
DEFAULT_OUTPUT = (
    PROJECT / "analysis" / "results" / "bgoegg_v1" /
    "bgo_threshold3_confirm_v1.json"
)


def train_and_apply(
    label: str,
    training_directory: Path,
    confirmation_directory: Path,
    feature_set: str,
    target: tuple[float, float],
) -> dict:
    features = HURDLE_FEATURE_SETS[feature_set]
    training = hurdle_feature_values(
        training_directory, label + ":training", features
    )
    model = fit_hurdle(training)
    validation_scores = {}
    for species, sample in training["samples"].items():
        mask = split_masks(sample["event_id"])["validation"]
        validation_scores[species] = hurdle_scores(model, sample, mask)
    cuts = choose_pairwise_thresholds(validation_scores, *target)
    validation = metrics_pairwise(validation_scores, cuts)

    confirmation = hurdle_feature_values(
        confirmation_directory, label + ":confirmation", features
    )
    confirmation_scores = {
        species: hurdle_scores(
            model, sample, np.ones(len(sample["event_id"]), dtype=bool)
        )
        for species, sample in confirmation["samples"].items()
    }
    return {
        "feature_set": feature_set,
        "classifier_inputs": features,
        "threshold_MeV": 3.0,
        "training_directory": str(training_directory.resolve()),
        "confirmation_directory": str(confirmation_directory.resolve()),
        "fit": "training eventID mod 4 == 0 only",
        "cut_selection": "training eventID mod 4 == 2 only",
        "target": {"pim_reject": target[0], "pi0_reject": target[1]},
        "threshold": [float(value) for value in cuts],
        "training_validation": validation,
        "confirmation_all_100k": metrics_pairwise(confirmation_scores, cuts),
        "confirmation_inputs": {
            species: sample["path"]
            for species, sample in confirmation["samples"].items()
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--confirm-cache", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    result = {
        "status": "complete",
        "analysis": "bgo_threshold3_confirm_v1",
        "confirmation_policy": (
            "feature set, fit sample, and validation cut are frozen from seed "
            "6302026; seed 7302026 is evaluated once without refitting"
        ),
        "ball_a20x40": train_and_apply(
            "ball_a20x40",
            args.cache / "normalized_standard" / "ball_a20x40",
            args.confirm_cache / "ball_a20x40",
            "energy_shape",
            TARGET,
        ),
        "egg29_zm10": train_and_apply(
            "egg29_zm10",
            args.cache / "standard" / "egg29zm10",
            args.confirm_cache / "egg29_zm10",
            "cluster_shape",
            GUARDED_TARGET,
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        key: value["confirmation_all_100k"]
        for key, value in result.items()
        if isinstance(value, dict) and "confirmation_all_100k" in value
    }, indent=2))
    print(args.output)


if __name__ == "__main__":
    main()
