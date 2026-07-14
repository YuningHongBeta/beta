#!/usr/bin/env python3
"""Explore BGO-only classifiers with a hard 3 MeV per-segment threshold.

All input BGO2 files must already have ``threshold_MeV == 3``.  The study
compares the existing full summary/pattern models with compact summary models
and an nHit-binned (hurdle) likelihood.  No energy below threshold is used.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from bgo_compare_v1 import (
    FEATURE_GROUPS,
    METHODS,
    analyze_geometry,
    choose_pairwise_thresholds,
    load_geometry,
    metrics_pairwise,
    split_masks,
)
from bgo_features_v2 import load_bgo2


PROJECT = Path(__file__).resolve().parent.parent
DEFAULT_CACHE = PROJECT / "tmp" / "bgo_threshold3_v1"
DEFAULT_OUTPUT = (
    PROJECT / "analysis" / "results" / "bgoegg_v1" /
    "bgo_threshold3_analysis_v1.json"
)
SPECIES_FILE = {"electron": "e.bgo2", "pim": "pim.bgo2", "pi0": "pi0.bgo2"}
TARGET = (0.92212, 0.98894)
# One-sided 95% binomial fluctuation guards for a 25k validation background
# sample, rounded upward.  This is an engineering operating point, not a new
# physics requirement.
GUARDED_TARGET = (0.9250, 0.9901)
THRESHOLD_MEV = 3.0
HIT_BINS = 6

# These small sets were declared after inspecting the T=3 validation behavior,
# but before the odd-sample values produced by this driver were inspected.
COMPACT_FEATURE_SETS = {
    "energy_rms": ["sumE", "allRmsDeg"],
    "energy_rms_isolation": ["sumE", "allRmsDeg", "isolatedEFrac4"],
    "energy_rms_hits": ["sumE", "allRmsDeg", "nHit"],
    "energy_rms_clusters": ["sumE", "allRmsDeg", "nCl4", "cl2OverCl1"],
    "energy_rms_local": ["sumE", "allRmsDeg", "local3x3Frac"],
    "energy_shape": [
        "sumE", "maxE", "allRmsDeg", "isolatedEFrac4", "local3x3Frac"
    ],
}

HURDLE_FEATURE_SETS = {
    "energy_rms": ["sumE", "allRmsDeg"],
    "energy_rms_isolation": ["sumE", "allRmsDeg", "isolatedEFrac4"],
    "energy_shape": [
        "sumE", "allRmsDeg", "isolatedEFrac4", "cl2Sum4",
        "local3x3Frac", "bgoLeadingTheta_deg",
    ],
    "cluster_shape": [
        "sumE", "allRmsDeg", "cl1MaxFrac4", "cl2OverCl1",
        "isolatedEFrac4", "local3x3Frac",
    ],
}

ADDITIONAL_STANDARD_GEOMETRIES = (
    "both3_28",
    "both3_28_zm10cm",
    "front3_back4_29",
    "ball_holes_29",
)


def triplet(directory: Path) -> dict[str, Path]:
    paths = {species: directory / filename for species, filename in SPECIES_FILE.items()}
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing BGO2 input: " + ", ".join(missing))
    return paths


def load_triplet(directory: Path, label: str, features: list[str] | None) -> dict:
    paths = triplet(directory)
    geometry = load_geometry(
        paths["electron"], paths["pim"], paths["pi0"], label, features
    )
    threshold = float(geometry["metadata"]["threshold_MeV"])
    if threshold != THRESHOLD_MEV:
        raise ValueError(f"{label}: expected {THRESHOLD_MEV} MeV, got {threshold}")
    return geometry


def selected_summary(result: dict) -> dict:
    target = next(iter(result["targets"].values()))
    return {
        "selected_method": target["selected_method"],
        "validation": target["methods"][target["selected_method"]]["validation"],
        "candidate_evaluation": target["selected_test"],
    }


def compact_scan(directory: Path, label: str) -> dict:
    candidates = {}
    for name, features in COMPACT_FEATURE_SETS.items():
        geometry = load_triplet(directory, f"{label}:compact:{name}", features)
        result = analyze_geometry(geometry, [TARGET], METHODS)
        candidates[name] = result
    selected_name = max(
        candidates,
        key=lambda name: selected_summary(candidates[name])["validation"]["e_keep"],
    )
    return {
        "selection_rule": "maximum validation electron retention across declared compact sets",
        "selected_feature_set": selected_name,
        "selected": selected_summary(candidates[selected_name]),
        "candidates": candidates,
    }


def hurdle_feature_values(directory: Path, label: str, features: list[str]) -> dict:
    # nHit controls the discrete component and is not repeated in the Gaussian.
    geometry = load_triplet(directory, label, ["nHit"] + features)
    for sample in geometry["samples"].values():
        names = sample["names"]
        sample["hit_bin"] = np.minimum(
            np.asarray(sample["raw"][:, names.index("nHit")], dtype=np.int64),
            HIT_BINS - 1,
        )
        values = np.asarray(
            sample["raw"][:, [names.index(name) for name in features]],
            dtype=np.float64,
        ).copy()
        # This reproduces the declared hurdle model: only total energy is
        # transformed; the remaining topology observables stay on their raw
        # reconstructed scale.
        if "sumE" in features:
            column = features.index("sumE")
            values[:, column] = np.log1p(np.maximum(values[:, column], 0.0))
        sample["conditional_x"] = values
    return geometry


def fit_hurdle(geometry: dict, pseudo_count: float = 50.0) -> dict:
    fit_parts = {}
    for species, sample in geometry["samples"].items():
        mask = split_masks(sample["event_id"])["fit"]
        fit_parts[species] = {
            "x": sample["conditional_x"][mask],
            "bin": sample["hit_bin"][mask],
        }
    pooled = np.vstack([part["x"] for part in fit_parts.values()])
    center = pooled.mean(axis=0)
    scale = pooled.std(axis=0)
    scale[scale < 1.0e-6] = 1.0
    species_models = {}
    for species, part in fit_parts.items():
        z = (part["x"] - center) / scale
        global_mu = z.mean(axis=0)
        global_var = z.var(axis=0)
        counts = np.bincount(part["bin"], minlength=HIT_BINS)
        bins = []
        for hit_bin in range(HIT_BINS):
            values = z[part["bin"] == hit_bin]
            nvalue = len(values)
            first = values.sum(axis=0) if nvalue else np.zeros(z.shape[1])
            second = (values * values).sum(axis=0) if nvalue else np.zeros(z.shape[1])
            denominator = nvalue + pseudo_count
            mean = (first + pseudo_count * global_mu) / denominator
            moment2 = (
                second + pseudo_count * (global_var + global_mu * global_mu)
            ) / denominator
            variance = np.maximum(moment2 - mean * mean, 0.08)
            bins.append({"mean": mean, "variance": variance})
        species_models[species] = {
            "log_prior": np.log((counts + 1.0) / (counts.sum() + HIT_BINS)),
            "counts": counts,
            "bins": bins,
        }
    return {
        "center": center,
        "scale": scale,
        "pseudo_count": pseudo_count,
        "species": species_models,
    }


def hurdle_log_likelihood(
    model: dict, species: str, values: np.ndarray, hit_bin: np.ndarray
) -> np.ndarray:
    z = (values - model["center"]) / model["scale"]
    species_model = model["species"][species]
    output = np.empty(len(values), dtype=np.float64)
    for bin_index in range(HIT_BINS):
        mask = hit_bin == bin_index
        if not np.any(mask):
            continue
        parameters = species_model["bins"][bin_index]
        delta = z[mask] - parameters["mean"]
        output[mask] = species_model["log_prior"][bin_index] - 0.5 * np.sum(
            delta * delta / parameters["variance"]
            + np.log(parameters["variance"]),
            axis=1,
        )
    return output


def hurdle_scores(model: dict, sample: dict, mask: np.ndarray) -> np.ndarray:
    values = sample["conditional_x"][mask]
    hit_bin = sample["hit_bin"][mask]
    log_e = hurdle_log_likelihood(model, "electron", values, hit_bin)
    return np.column_stack(
        [
            hurdle_log_likelihood(model, "pim", values, hit_bin) - log_e,
            hurdle_log_likelihood(model, "pi0", values, hit_bin) - log_e,
        ]
    )


def acceptance_by_hit_bin(
    scores: dict[str, np.ndarray], bins: dict[str, np.ndarray], cuts: tuple[float, float]
) -> dict:
    output = {}
    for species, values in scores.items():
        accepted = (values[:, 0] < cuts[0]) & (values[:, 1] < cuts[1])
        rows = {}
        for hit_bin in range(HIT_BINS):
            mask = bins[species] == hit_bin
            rows[str(hit_bin) if hit_bin < HIT_BINS - 1 else "5+"] = {
                "count": int(mask.sum()),
                "acceptance": float(accepted[mask].mean()) if np.any(mask) else None,
            }
        output[species] = rows
    return output


def analyze_hurdle(
    directory: Path,
    label: str,
    features: list[str],
    target: tuple[float, float] = TARGET,
) -> dict:
    geometry = hurdle_feature_values(directory, label, features)
    model = fit_hurdle(geometry)
    masks = {
        species: split_masks(sample["event_id"])
        for species, sample in geometry["samples"].items()
    }
    scores = {
        split: {
            species: hurdle_scores(model, sample, masks[species][split])
            for species, sample in geometry["samples"].items()
        }
        for split in ("validation", "test")
    }
    cuts = choose_pairwise_thresholds(scores["validation"], *target)
    test_bins = {
        species: sample["hit_bin"][masks[species]["test"]]
        for species, sample in geometry["samples"].items()
    }
    return {
        "classifier_inputs": features,
        "discrete_input": "nHit capped into 0,1,2,3,4,5+ bins",
        "conditional_model": "diagonal Gaussian shrunk to species-wide moments",
        "pseudo_count": model["pseudo_count"],
        "target": {"pim_reject": target[0], "pi0_reject": target[1]},
        "fit_bin_counts": {
            species: [int(value) for value in parameters["counts"]]
            for species, parameters in model["species"].items()
        },
        "threshold": [float(value) for value in cuts],
        "validation": metrics_pairwise(scores["validation"], cuts),
        "candidate_evaluation": metrics_pairwise(scores["test"], cuts),
        "test_acceptance_by_hit_bin": acceptance_by_hit_bin(
            scores["test"], test_bins, cuts
        ),
    }


def hurdle_scan(directory: Path, label: str) -> dict:
    candidates = {
        name: analyze_hurdle(directory, f"{label}:hurdle:{name}", features)
        for name, features in HURDLE_FEATURE_SETS.items()
    }
    selected_name = max(
        candidates,
        key=lambda name: candidates[name]["validation"]["e_keep"],
    )
    guarded = analyze_hurdle(
        directory,
        f"{label}:hurdle:{selected_name}:guarded",
        HURDLE_FEATURE_SETS[selected_name],
        GUARDED_TARGET,
    )
    return {
        "selection_rule": "maximum validation electron retention across declared hurdle sets",
        "selected_feature_set": selected_name,
        "selected": candidates[selected_name],
        "guarded_operating_point": guarded,
        "candidates": candidates,
    }


def ncluster_distribution(directory: Path, label: str) -> dict:
    geometry = load_triplet(directory, label, ["nCl4"])
    output = {}
    for species, sample in geometry["samples"].items():
        values = np.asarray(sample["x"][:, 0], dtype=np.int64)
        output[species] = {
            "count": int(len(values)),
            "ncluster_0": float(np.mean(values == 0)),
            "ncluster_1": float(np.mean(values == 1)),
            "ncluster_ge_2": float(np.mean(values >= 2)),
        }
    return output


def two_species_ncluster(directory: Path, threshold: int) -> dict:
    output = {}
    for species, filename in (("electron", "e"), ("pim", "pim")):
        path = directory / f"{filename}_T{threshold}.bgo2"
        data, meta, manifest = load_bgo2(path)
        if float(meta["threshold_MeV"]) != float(threshold):
            raise ValueError(f"{path}: threshold metadata mismatch")
        names = list(manifest["features"])
        values = np.asarray(data[:, names.index("nCl4")], dtype=np.int64)
        output[species] = {
            "count": int(len(values)),
            "ncluster_0": float(np.mean(values == 0)),
            "ncluster_1": float(np.mean(values == 1)),
            "ncluster_ge_2": float(np.mean(values >= 2)),
        }
    return output


def hong_thesis_check(directory: Path) -> dict:
    # Hong thesis Table 6.7, uniform-theta 15x15, Geant4 11.2.2 INCLXX MOD.
    thesis = {
        "T1": {
            "electron": {"ncluster_0": 0.008, "ncluster_1": 0.757,
                         "ncluster_ge_2": 0.235},
            "pim": {"ncluster_0": 0.011, "ncluster_1": 0.118,
                    "ncluster_ge_2": 0.872},
        },
        "T3": {
            "electron": {"ncluster_0": 0.011, "ncluster_1": 0.904,
                         "ncluster_ge_2": 0.085},
            "pim": {"ncluster_0": 0.017, "ncluster_1": 0.274,
                    "ncluster_ge_2": 0.709},
        },
    }
    current = {tag: two_species_ncluster(directory, threshold)
               for tag, threshold in (("T1", 1), ("T3", 3))}
    differences = {}
    for tag in thesis:
        differences[tag] = {}
        for species in thesis[tag]:
            differences[tag][species] = {
                key + "_percentage_point":
                    100.0 * (current[tag][species][key] - value)
                for key, value in thesis[tag][species].items()
            }
    return {
        "reference": (
            "260401_MThesis_YuningHong.pdf Table 6.7; segment-threshold "
            "columns, uniform-theta 15x15, four-neighbor, N=100000"
        ),
        "comparability": (
            "same segmentation/threshold/neighbor/count; current ROOT is a "
            "newer independent production and is not the thesis event sample"
        ),
        "thesis": thesis,
        "current": current,
        "current_minus_thesis": differences,
    }


def analyze_one(
    label: str, standard_directory: Path, pattern_directory: Path
) -> dict:
    standard = analyze_geometry(
        load_triplet(standard_directory, label + ":standard", FEATURE_GROUPS["bgo"]),
        [TARGET],
        METHODS,
    )
    pattern = analyze_geometry(
        load_triplet(pattern_directory, label + ":pattern", None),
        [TARGET],
        METHODS[:-1],
    )
    return {
        "standard_full": standard,
        "standard_full_selected": selected_summary(standard),
        "pattern_full": pattern,
        "pattern_full_selected": selected_summary(pattern),
        "compact": compact_scan(standard_directory, label),
        "hurdle": hurdle_scan(standard_directory, label),
        "ncluster_full_sample": ncluster_distribution(
            standard_directory, label + ":ncluster"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    parser.add_argument(
        "--hong-cache", type=Path,
        default=PROJECT / "tmp" / "hong_thesis_check" / "u15x15",
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    standard = {
        "ball_a20x40": args.cache.parent / "bgoegg_study_v1" / "standard" /
        "ball__a20x40",
        "egg22": args.cache.parent / "bgoegg_study_v1" / "standard" /
        "production__egg_none",
        "egg29_zm10": args.cache / "standard" / "egg29zm10",
    }
    # The older cache names include the threshold in the filename.  Build a
    # normalized read-only view by accepting either naming convention.
    normalized = args.cache / "normalized_standard"
    for label in ("ball_a20x40", "egg22"):
        destination = normalized / label
        destination.mkdir(parents=True, exist_ok=True)
        for species, filename in SPECIES_FILE.items():
            source_name = filename.replace(".bgo2", "_T3.bgo2")
            source = standard[label] / source_name
            target = destination / filename
            side_source = Path(str(source) + ".json")
            side_target = Path(str(target) + ".json")
            if target.exists() or target.is_symlink():
                target.unlink()
            if side_target.exists() or side_target.is_symlink():
                side_target.unlink()
            target.symlink_to(source.resolve())
            side_target.symlink_to(side_source.resolve())
        standard[label] = destination
    pattern = {
        "ball_a20x40": args.cache / "pattern" / "ball",
        "egg22": args.cache / "pattern" / "egg22",
        "egg29_zm10": args.cache / "pattern" / "egg29zm10",
    }
    geometries = {
        label: analyze_one(label, standard[label], pattern[label])
        for label in standard
    }
    additional = {}
    for label in ADDITIONAL_STANDARD_GEOMETRIES:
        directory = args.cache / "standard_scan" / label
        additional[label] = {
            "compact": compact_scan(directory, label),
            "hurdle": hurdle_scan(directory, label),
            "ncluster_full_sample": ncluster_distribution(
                directory, label + ":ncluster"
            ),
        }
    result = {
        "status": "complete",
        "analysis": "bgo_threshold3_analysis_v1",
        "threshold_MeV": THRESHOLD_MEV,
        "threshold_policy": "hard per-segment cut; deposits below 3 MeV are not classifier inputs",
        "target": {"pim_reject": TARGET[0], "pi0_reject": TARGET[1]},
        "split": {
            "fit": "eventID mod 4 == 0",
            "validation": "eventID mod 4 == 2; method/feature/cut selection",
            "candidate_evaluation": "odd eventID; not a study-wide blind holdout",
        },
        "geometries": geometries,
        "additional_standard_geometry_scan": additional,
        "hong_thesis_check": hong_thesis_check(args.hong_cache),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    concise = {}
    for label, geometry in geometries.items():
        concise[label] = {
            "standard": geometry["standard_full_selected"],
            "pattern": geometry["pattern_full_selected"],
            "compact": {
                "feature_set": geometry["compact"]["selected_feature_set"],
                **geometry["compact"]["selected"],
            },
            "hurdle": {
                "feature_set": geometry["hurdle"]["selected_feature_set"],
                **geometry["hurdle"]["selected"]["candidate_evaluation"],
            },
            "hurdle_guarded": {
                "feature_set": geometry["hurdle"]["selected_feature_set"],
                **geometry["hurdle"]["guarded_operating_point"][
                    "candidate_evaluation"
                ],
            },
        }
    print(json.dumps(concise, indent=2))
    print(args.output)


if __name__ == "__main__":
    main()
