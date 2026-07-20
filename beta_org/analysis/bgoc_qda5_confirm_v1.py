#!/usr/bin/env python3
"""Freeze and confirm the five-variable BGOC T3 QDA on an independent seed."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np

from bgo_compare_v1 import (
    choose_threshold,
    fit_method,
    load_geometry,
    metrics,
    score_method,
    split_masks,
)
from bgo_optimize_v2 import BGO_FEATURES


PROJECT = Path(__file__).resolve().parent.parent
FEATURES = [
    "sumE",
    "maxE",
    "allRmsDeg",
    "isolatedEFrac4",
    "local3x3Frac",
]
SPECIES_FILE = {"electron": "e.bgo2", "pim": "pim.bgo2", "pi0": "pi0.bgo2"}
PRIMARY = {"electron": "e", "pim": "pim", "pi0": "pi0"}
TARGET = (0.88448, 0.98448)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def triplet(directory: Path) -> dict[str, Path]:
    paths = {species: directory / name for species, name in SPECIES_FILE.items()}
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing BGO2 input: " + ", ".join(missing))
    return paths


def geometry(directory: Path, label: str, features: list[str] | None = None) -> dict:
    paths = triplet(directory)
    return load_geometry(
        paths["electron"], paths["pim"], paths["pi0"], label, features or FEATURES
    )


def extract_confirmation(args: argparse.Namespace) -> dict[str, Path]:
    args.cache.mkdir(parents=True, exist_ok=True)
    roots = {}
    for species, primary in PRIMARY.items():
        root = args.root_directory / f"{args.tag}_u15x15_{primary}.root"
        if not root.is_file():
            raise FileNotFoundError(root)
        output = args.cache / SPECIES_FILE[species]
        command = [str(args.extractor), str(root), str(output), "3"]
        subprocess.run(command, check=True)
        roots[species] = root
    return roots


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--training",
        type=Path,
        default=PROJECT / "tmp/bgoc_u15x15_feasibility_v1/bgo2/none/training_T3",
    )
    parser.add_argument(
        "--root-directory", type=Path, default=PROJECT / "build_opt/output/scan"
    )
    parser.add_argument(
        "--tag", default="bgoc_u15x15_qda5_confirm_s9302026"
    )
    parser.add_argument(
        "--extractor", type=Path, default=PROJECT / "analysis/bgo_extract_features_v2"
    )
    parser.add_argument(
        "--cache", type=Path, default=PROJECT / "tmp/bgoc_qda5_confirm_s9302026"
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=PROJECT / "analysis/results/bgoegg_v1/bgoc_qda5_confirm_s9302026.json",
    )
    parser.add_argument("--skip-extract", action="store_true")
    args = parser.parse_args()

    roots = {}
    if args.skip_extract:
        triplet(args.cache)
        for species, primary in PRIMARY.items():
            roots[species] = args.root_directory / f"{args.tag}_u15x15_{primary}.root"
    else:
        roots = extract_confirmation(args)

    training = geometry(args.training, "seed6302026-training-T3")
    confirmation = geometry(args.cache, "seed9302026-confirmation-T3")
    for loaded in (training, confirmation):
        if float(loaded["metadata"]["threshold_MeV"]) != 3.0:
            raise ValueError("all inputs must use a 3 MeV hard per-cell threshold")

    fit_parts = {
        species: sample["x"][split_masks(sample["event_id"])["fit"]]
        for species, sample in training["samples"].items()
    }
    model = fit_method("qda", fit_parts)
    validation_scores = {
        species: score_method(
            "qda", model, sample["x"][split_masks(sample["event_id"])["validation"]]
        )
        for species, sample in training["samples"].items()
    }
    threshold = choose_threshold(validation_scores, *TARGET)
    validation = metrics(validation_scores, threshold)
    confirmation_scores = {
        species: score_method("qda", model, sample["x"])
        for species, sample in confirmation["samples"].items()
    }
    confirmation_result = metrics(confirmation_scores, threshold)

    expected_validation = {
        "e_keep": 0.60148,
        "pim_reject": 0.88452,
        "pi0_reject": 0.98716,
    }
    for name, expected in expected_validation.items():
        if validation[name] != expected:
            raise RuntimeError(
                f"frozen validation mismatch for {name}: {validation[name]} != {expected}"
            )

    reference_training = geometry(
        args.training, "seed6302026-training-T3-quadratic-reference", BGO_FEATURES
    )
    reference_confirmation = geometry(
        args.cache, "seed9302026-confirmation-T3-quadratic-reference", BGO_FEATURES
    )
    reference_fit = {
        species: sample["x"][split_masks(sample["event_id"])["fit"]]
        for species, sample in reference_training["samples"].items()
    }
    reference_model = fit_method("quadratic_ridge", reference_fit)
    reference_validation_scores = {
        species: score_method(
            "quadratic_ridge",
            reference_model,
            sample["x"][split_masks(sample["event_id"])["validation"]],
        )
        for species, sample in reference_training["samples"].items()
    }
    reference_threshold = choose_threshold(reference_validation_scores, *TARGET)
    reference_validation = metrics(reference_validation_scores, reference_threshold)
    reference_confirmation_scores = {
        species: score_method("quadratic_ridge", reference_model, sample["x"])
        for species, sample in reference_confirmation["samples"].items()
    }
    reference_confirmation_result = metrics(
        reference_confirmation_scores, reference_threshold
    )
    expected_reference_validation = {
        "e_keep": 0.61076,
        "pim_reject": 0.88452,
        "pi0_reject": 0.98688,
    }
    for name, expected in expected_reference_validation.items():
        if reference_validation[name] != expected:
            raise RuntimeError(
                "frozen quadratic-ridge validation mismatch for "
                f"{name}: {reference_validation[name]} != {expected}"
            )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    result = {
        "status": "complete",
        "analysis": "bgoc_qda5_confirm_v1",
        "decision": "five-variable QDA frozen after the seed7302026 audit",
        "geometry": "current BGOC 15x15 uniform-theta, z=0",
        "bgo_cell_threshold_MeV": 3.0,
        "model": {
            "method": "QDA",
            "features": FEATURES,
            "energy_transform": "log1p for sumE and maxE",
            "standardization": "pooled fit-sample mean and standard deviation",
            "covariance_regularization": 0.20,
            "score": "log((L_pim + L_pi0)/2) - log(L_e)",
            "acceptance": "score < threshold",
        },
        "split": {
            "fit": "seed6302026 eventID mod 4 == 0; 25,000/species",
            "cut_selection": "seed6302026 eventID mod 4 == 2; 25,000/species",
            "confirmation": "seed9302026; all 100,000/species",
        },
        "target": {"pim_reject": TARGET[0], "pi0_reject": TARGET[1]},
        "threshold": float(threshold),
        "training_validation": validation,
        "confirmation_all": confirmation_result,
        "reference_quadratic_ridge": {
            "role": "previous 23-variable model evaluated on the same new seed; not used to retune QDA5",
            "features": BGO_FEATURES,
            "threshold": float(reference_threshold),
            "training_validation": reference_validation,
            "confirmation_all": reference_confirmation_result,
        },
        "source": {
            "training_bgo2": {
                species: {"path": str(path), "sha256": sha256(path)}
                for species, path in triplet(args.training).items()
            },
            "confirmation_root": {
                species: {"path": str(path), "sha256": sha256(path)}
                for species, path in roots.items()
            },
            "confirmation_bgo2": {
                species: {"path": str(path), "sha256": sha256(path)}
                for species, path in triplet(args.cache).items()
            },
            "extractor": {"path": str(args.extractor), "sha256": sha256(args.extractor)},
        },
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(args.output)
    print(json.dumps(confirmation_result, indent=2))


if __name__ == "__main__":
    main()
