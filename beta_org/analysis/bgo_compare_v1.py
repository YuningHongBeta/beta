#!/usr/bin/env python3
"""Compare BGO ball and BGOegg analysis performance on independent test data.

The input files are compact BGO2 feature files made by
``bgo_extract_features_v2``.  All classifiers use reconstructed BGO shower
features only.  Event IDs divisible by four fit the score model, event IDs
equal to two modulo four set the operating threshold and choose the method,
and odd event IDs are touched only for the final test measurement.
"""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path

for _variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_variable] = "1"

import numpy as np

from bgo_features_v2 import load_bgo2
from bgo_optimize_v2 import BGO_FEATURES, LOG_FEATURES


METHODS = ("fisher", "qda", "histogram", "pairwise_qda", "quadratic_ridge")
FEATURE_GROUPS = {
    "bgo": BGO_FEATURES,
    "bgo_tlc": BGO_FEATURES + ["tlcAnySegHitGe1_eff05"],
    "bgo_th_tlc": BGO_FEATURES
    + [
        "thMatchedDEdx_MeV_per_mm",
        "absDeltaZ_mm",
        "absDeltaZValid",
        "absDeltaZLt90",
        "tlcAnySegHitGe1_eff05",
    ],
    "all": None,
}


def wilson(value: float, n: int) -> list[float]:
    z = 1.959963984540054
    denominator = 1.0 + z * z / n
    center = (value + z * z / (2.0 * n)) / denominator
    half = z * math.sqrt(
        value * (1.0 - value) / n + z * z / (4.0 * n * n)
    ) / denominator
    return [float(center - half), float(center + half)]


def load_species(path: Path, feature_names: list[str] | None) -> dict:
    data, meta, manifest = load_bgo2(path)
    names = list(manifest["features"])
    if feature_names is None:
        feature_names = [name for name in names if name != "eventID"]
    missing = sorted(set(feature_names + ["eventID"]) - set(names))
    if missing:
        raise ValueError(f"{path}: missing features {missing}")
    indices = [names.index(name) for name in feature_names]
    values = np.asarray(data[:, indices], dtype=np.float64)
    for column, name in enumerate(feature_names):
        if (
            name in LOG_FEATURES
            or "SumE_" in name
            or "MaxE_" in name
            or name.startswith("patchE_")
        ):
            values[:, column] = np.log1p(np.maximum(0.0, values[:, column]))
    event_id = np.asarray(data[:, names.index("eventID")], dtype=np.int64)
    if np.any(event_id != data[:, names.index("eventID")]):
        raise ValueError(f"{path}: non-integral eventID")
    return {
        "path": str(path),
        "x": values,
        "event_id": event_id,
        "meta": meta,
        "manifest": manifest,
        "raw": data,
        "names": names,
        "classifier_inputs": feature_names,
    }


def comparable_meta(sample: dict) -> dict:
    meta = sample["meta"]
    comparable = {
        key: meta[key]
        for key in (
            "version", "ncol", "nLayer", "nSector", "segmentationMode",
            "physicsFlag", "thetaMin_deg", "thetaMax_deg", "threshold_MeV",
        )
    }
    for key in ("geometryMode", "rMin_cm", "thickness_cm", "bgoZOffset_cm"):
        comparable[key] = sample["manifest"].get(key)
    return comparable


def load_geometry(
    electron: Path,
    pion: Path,
    pi0: Path,
    label: str,
    feature_names: list[str] | None,
) -> dict:
    electron_sample = load_species(electron, feature_names)
    selected_features = electron_sample["classifier_inputs"]
    samples = {
        "electron": electron_sample,
        "pim": load_species(pion, selected_features),
        "pi0": load_species(pi0, selected_features),
    }
    reference = comparable_meta(samples["electron"])
    for species, sample in samples.items():
        if comparable_meta(sample) != reference:
            raise ValueError(
                f"{label}: metadata mismatch for {species}: "
                f"{comparable_meta(sample)} != {reference}"
            )
        if len(np.unique(sample["event_id"])) != len(sample["event_id"]):
            raise ValueError(f"{label}: duplicate {species} eventID")
    return {
        "label": label,
        "samples": samples,
        "metadata": reference,
        "feature_names": selected_features,
    }


def split_masks(event_id: np.ndarray) -> dict[str, np.ndarray]:
    return {
        "fit": (event_id & 3) == 0,
        "validation": (event_id & 3) == 2,
        "test": (event_id & 1) == 1,
    }


def standardizer(parts: list[np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    pooled = np.vstack(parts)
    center = pooled.mean(axis=0)
    scale = pooled.std(axis=0)
    scale[scale < 1.0e-10] = 1.0
    return center, scale


def fit_fisher(parts: dict[str, np.ndarray], regularization: float = 0.08) -> dict:
    center, scale = standardizer(list(parts.values()))
    z = {name: (values - center) / scale for name, values in parts.items()}
    background = np.vstack([z["pim"], z["pi0"]])
    ce = np.atleast_2d(np.cov(z["electron"], rowvar=False))
    cb = np.atleast_2d(np.cov(background, rowvar=False))
    covariance = 0.5 * (ce + cb)
    covariance = (
        (1.0 - regularization) * covariance
        + regularization * np.eye(covariance.shape[0])
    )
    weights = np.linalg.solve(
        covariance, background.mean(axis=0) - z["electron"].mean(axis=0)
    )
    return {
        "center": center,
        "scale": scale,
        "weights": weights,
        "regularization": regularization,
    }


def score_fisher(model: dict, values: np.ndarray) -> np.ndarray:
    return ((values - model["center"]) / model["scale"]) @ model["weights"]


def regularized_gaussian(values: np.ndarray, regularization: float) -> dict:
    mean = values.mean(axis=0)
    covariance = np.atleast_2d(np.cov(values, rowvar=False))
    diagonal = np.diag(np.diag(covariance))
    covariance = (1.0 - regularization) * covariance + regularization * diagonal
    floor = max(float(np.trace(covariance) / len(mean)) * 1.0e-6, 1.0e-8)
    covariance += floor * np.eye(len(mean))
    sign, logdet = np.linalg.slogdet(covariance)
    if sign <= 0:
        raise RuntimeError("QDA covariance is not positive definite")
    return {
        "mean": mean,
        "precision": np.linalg.inv(covariance),
        "logdet": float(logdet),
    }


def gaussian_log_density(model: dict, values: np.ndarray) -> np.ndarray:
    delta = values - model["mean"]
    return -0.5 * (
        np.einsum("ij,jk,ik->i", delta, model["precision"], delta)
        + model["logdet"]
    )


def fit_qda(parts: dict[str, np.ndarray], regularization: float = 0.20) -> dict:
    center, scale = standardizer(list(parts.values()))
    z = {name: (values - center) / scale for name, values in parts.items()}
    return {
        "center": center,
        "scale": scale,
        "species": {
            name: regularized_gaussian(values, regularization)
            for name, values in z.items()
        },
        "regularization": regularization,
    }


def score_qda(model: dict, values: np.ndarray) -> np.ndarray:
    z = (values - model["center"]) / model["scale"]
    log_e = gaussian_log_density(model["species"]["electron"], z)
    log_pim = gaussian_log_density(model["species"]["pim"], z)
    log_pi0 = gaussian_log_density(model["species"]["pi0"], z)
    log_background = np.logaddexp(log_pim, log_pi0) - math.log(2.0)
    return log_background - log_e


def score_pairwise_qda(model: dict, values: np.ndarray) -> np.ndarray:
    """Return separate pi-minus/electron and pi-zero/electron log odds."""
    z = (values - model["center"]) / model["scale"]
    log_e = gaussian_log_density(model["species"]["electron"], z)
    return np.column_stack(
        [
            gaussian_log_density(model["species"]["pim"], z) - log_e,
            gaussian_log_density(model["species"]["pi0"], z) - log_e,
        ]
    )


def fit_histogram(
    parts: dict[str, np.ndarray], nbin: int = 32, smoothing: float = 1.0
) -> dict:
    pooled = np.vstack(list(parts.values()))
    quantiles = np.linspace(0.0, 1.0, nbin + 1)
    edges = []
    log_probabilities = {name: [] for name in parts}
    for column in range(pooled.shape[1]):
        edge = np.unique(np.quantile(pooled[:, column], quantiles))
        if len(edge) < 3:
            edge = np.array([-np.inf, np.inf])
        else:
            edge[0], edge[-1] = -np.inf, np.inf
        edges.append(edge)
        nb = len(edge) - 1
        for name, values in parts.items():
            index = np.searchsorted(edge, values[:, column], side="right") - 1
            index = np.clip(index, 0, nb - 1)
            counts = np.bincount(index, minlength=nb).astype(np.float64)
            probabilities = (counts + smoothing) / (counts.sum() + smoothing * nb)
            log_probabilities[name].append(np.log(probabilities))
    return {
        "edges": edges,
        "log_probabilities": log_probabilities,
        "nbin_requested": nbin,
        "smoothing": smoothing,
    }


def histogram_log_density(model: dict, species: str, values: np.ndarray) -> np.ndarray:
    result = np.zeros(len(values), dtype=np.float64)
    for column, edge in enumerate(model["edges"]):
        nb = len(edge) - 1
        index = np.searchsorted(edge, values[:, column], side="right") - 1
        index = np.clip(index, 0, nb - 1)
        result += model["log_probabilities"][species][column][index]
    return result


def score_histogram(model: dict, values: np.ndarray) -> np.ndarray:
    log_e = histogram_log_density(model, "electron", values)
    log_pim = histogram_log_density(model, "pim", values)
    log_pi0 = histogram_log_density(model, "pi0", values)
    return np.logaddexp(log_pim, log_pi0) - math.log(2.0) - log_e


def quadratic_features(z: np.ndarray) -> np.ndarray:
    nrow, ncol = z.shape
    nout = 1 + ncol + ncol * (ncol + 1) // 2
    output = np.empty((nrow, nout), dtype=np.float64)
    output[:, 0] = 1.0
    output[:, 1 : 1 + ncol] = z
    offset = 1 + ncol
    for first in range(ncol):
        width = ncol - first
        output[:, offset : offset + width] = z[:, first, None] * z[:, first:]
        offset += width
    return output


def fit_quadratic_ridge(
    parts: dict[str, np.ndarray], regularization: float = 0.25
) -> dict:
    center, scale = standardizer(list(parts.values()))
    design_parts = []
    target_parts = []
    weight_parts = []
    for species, values in parts.items():
        design = quadratic_features((values - center) / scale)
        design_parts.append(design)
        target_parts.append(
            np.full(len(values), 0.0 if species == "electron" else 1.0)
        )
        # Equal total weight for electron, pi-minus, and pi-zero.
        weight_parts.append(np.full(len(values), 1.0 / len(values)))
    design = np.vstack(design_parts)
    target = np.concatenate(target_parts)
    weight = np.concatenate(weight_parts)
    xtwx = design.T @ (design * weight[:, None])
    penalty = regularization * np.eye(xtwx.shape[0])
    penalty[0, 0] = 0.0
    coefficients = np.linalg.solve(
        xtwx + penalty, design.T @ (weight * target)
    )
    return {
        "center": center,
        "scale": scale,
        "coefficients": coefficients,
        "regularization": regularization,
    }


def score_quadratic_ridge(model: dict, values: np.ndarray) -> np.ndarray:
    z = (values - model["center"]) / model["scale"]
    return quadratic_features(z) @ model["coefficients"]


def fit_method(method: str, parts: dict[str, np.ndarray]) -> dict:
    if method == "fisher":
        return fit_fisher(parts)
    if method == "qda":
        return fit_qda(parts)
    if method == "pairwise_qda":
        return fit_qda(parts)
    if method == "histogram":
        return fit_histogram(parts)
    if method == "quadratic_ridge":
        return fit_quadratic_ridge(parts)
    raise ValueError(method)


def score_method(method: str, model: dict, values: np.ndarray) -> np.ndarray:
    if method == "fisher":
        return score_fisher(model, values)
    if method == "qda":
        return score_qda(model, values)
    if method == "pairwise_qda":
        return score_pairwise_qda(model, values)
    if method == "histogram":
        return score_histogram(model, values)
    if method == "quadratic_ridge":
        return score_quadratic_ridge(model, values)
    raise ValueError(method)


def rejection_threshold(scores: np.ndarray, target: float) -> float:
    # Reject score >= threshold.  The selected order statistic guarantees at
    # least the requested rejection on the threshold-setting sample.
    index = max(0, int(math.floor((1.0 - target) * len(scores))))
    return float(np.partition(scores, index)[index])


def choose_threshold(
    validation_scores: dict[str, np.ndarray], target_pim: float, target_pi0: float
) -> float:
    return min(
        rejection_threshold(validation_scores["pim"], target_pim),
        rejection_threshold(validation_scores["pi0"], target_pi0),
    )


def metrics(scores: dict[str, np.ndarray], threshold: float) -> dict:
    e_keep = float(np.mean(scores["electron"] < threshold))
    pim_reject = float(np.mean(scores["pim"] >= threshold))
    pi0_reject = float(np.mean(scores["pi0"] >= threshold))
    return {
        "e_keep": e_keep,
        "pim_reject": pim_reject,
        "pi0_reject": pi0_reject,
        "e_keep_95ci": wilson(e_keep, len(scores["electron"])),
        "pim_reject_95ci": wilson(pim_reject, len(scores["pim"])),
        "pi0_reject_95ci": wilson(pi0_reject, len(scores["pi0"])),
        "counts": {name: int(len(value)) for name, value in scores.items()},
    }


def choose_pairwise_thresholds(
    validation_scores: dict[str, np.ndarray],
    target_pim: float,
    target_pi0: float,
    nscan: int = 512,
) -> tuple[float, float]:
    """Maximize validation electron acceptance in a two-score rectangle.

    An event is accepted only when both pairwise background/electron log odds
    are below their thresholds.  The scan uses only validation events.  For
    each first-score threshold, the largest allowed second-score threshold is
    obtained from the two background order statistics.
    """
    pooled_first = np.concatenate(
        [values[:, 0] for values in validation_scores.values()]
    )
    quantiles = np.linspace(0.0, 1.0, nscan + 1)
    candidates = np.unique(np.quantile(pooled_first, quantiles))
    candidates = np.concatenate(([-np.inf], candidates, [np.inf]))
    allowed = {
        "pim": int(math.floor((1.0 - target_pim) * len(validation_scores["pim"]))),
        "pi0": int(math.floor((1.0 - target_pi0) * len(validation_scores["pi0"]))),
    }
    best = (-1.0, -np.inf, -np.inf)
    electron = validation_scores["electron"]
    for threshold_first in candidates:
        second_limits = []
        for species in ("pim", "pi0"):
            values = validation_scores[species]
            subset = values[values[:, 0] < threshold_first, 1]
            if len(subset) <= allowed[species]:
                second_limits.append(np.inf)
            else:
                second_limits.append(
                    float(np.partition(subset, allowed[species])[allowed[species]])
                )
        threshold_second = min(second_limits)
        keep = float(
            np.mean(
                (electron[:, 0] < threshold_first)
                & (electron[:, 1] < threshold_second)
            )
        )
        candidate = (keep, float(threshold_first), float(threshold_second))
        if candidate > best:
            best = candidate
    return best[1], best[2]


def metrics_pairwise(
    scores: dict[str, np.ndarray], thresholds: tuple[float, float]
) -> dict:
    accepted = {
        species: (values[:, 0] < thresholds[0]) & (values[:, 1] < thresholds[1])
        for species, values in scores.items()
    }
    e_keep = float(accepted["electron"].mean())
    pim_reject = float(1.0 - accepted["pim"].mean())
    pi0_reject = float(1.0 - accepted["pi0"].mean())
    return {
        "e_keep": e_keep,
        "pim_reject": pim_reject,
        "pi0_reject": pi0_reject,
        "e_keep_95ci": wilson(e_keep, len(accepted["electron"])),
        "pim_reject_95ci": wilson(pim_reject, len(accepted["pim"])),
        "pi0_reject_95ci": wilson(pi0_reject, len(accepted["pi0"])),
        "counts": {name: int(len(value)) for name, value in scores.items()},
    }


def cluster_baseline(geometry: dict) -> dict:
    if any("nCl4" not in sample["names"] for sample in geometry["samples"].values()):
        return None
    output = {}
    for split in ("validation", "test"):
        values = {}
        for species, sample in geometry["samples"].items():
            mask = split_masks(sample["event_id"])[split]
            ncluster = np.asarray(
                sample["raw"][mask, sample["names"].index("nCl4")]
            )
            values[species] = ncluster == 1
        e_keep = float(values["electron"].mean())
        pim_reject = float(1.0 - values["pim"].mean())
        pi0_reject = float(1.0 - values["pi0"].mean())
        output[split] = {
            "e_keep": e_keep,
            "pim_reject": pim_reject,
            "pi0_reject": pi0_reject,
        }
        if split == "test":
            output[split].update(
                {
                    "e_keep_95ci": wilson(e_keep, len(values["electron"])),
                    "pim_reject_95ci": wilson(
                        pim_reject, len(values["pim"])
                    ),
                    "pi0_reject_95ci": wilson(pi0_reject, len(values["pi0"])),
                }
            )
    return output


def analyze_geometry(
    geometry: dict, targets: list[tuple[float, float]], methods: tuple[str, ...]
) -> dict:
    masks = {
        species: split_masks(sample["event_id"])
        for species, sample in geometry["samples"].items()
    }
    fit_parts = {
        species: sample["x"][masks[species]["fit"]]
        for species, sample in geometry["samples"].items()
    }
    models = {}
    for method in methods:
        if method == "pairwise_qda" and "qda" in models:
            models[method] = models["qda"]
        else:
            models[method] = fit_method(method, fit_parts)
    scores = {
        method: {
            split: {
                species: score_method(
                    method, model, sample["x"][masks[species][split]]
                )
                for species, sample in geometry["samples"].items()
            }
            for split in ("validation", "test")
        }
        for method, model in models.items()
    }
    target_results = {}
    for target_pim, target_pi0 in targets:
        key = f"pim{100*target_pim:.1f}_pi0{100*target_pi0:.1f}"
        method_results = {}
        for method in methods:
            if method == "pairwise_qda":
                threshold = choose_pairwise_thresholds(
                    scores[method]["validation"], target_pim, target_pi0
                )
                validation_metrics = metrics_pairwise(
                    scores[method]["validation"], threshold
                )
                test_metrics = metrics_pairwise(scores[method]["test"], threshold)
            else:
                threshold = choose_threshold(
                    scores[method]["validation"], target_pim, target_pi0
                )
                validation_metrics = metrics(
                    scores[method]["validation"], threshold
                )
                test_metrics = metrics(scores[method]["test"], threshold)
            method_results[method] = {
                "threshold": threshold,
                "validation": validation_metrics,
                "test": test_metrics,
            }
        selected = max(
            methods,
            key=lambda method: method_results[method]["validation"]["e_keep"],
        )
        target_results[key] = {
            "target_pim_reject": target_pim,
            "target_pi0_reject": target_pi0,
            "selection_rule": "maximum validation e_keep; test not inspected",
            "selected_method": selected,
            "methods": method_results,
            "selected_test": method_results[selected]["test"],
        }
    return {
        "label": geometry["label"],
        "classifier_inputs": geometry["feature_names"],
        "metadata": geometry["metadata"],
        "source_bgo2": {
            species: sample["path"] for species, sample in geometry["samples"].items()
        },
        "split_counts": {
            species: {
                split: int(mask.sum())
                for split, mask in masks[species].items()
            }
            for species in masks
        },
        "cluster_ncl4_eq_1": cluster_baseline(geometry),
        "targets": target_results,
    }


def parse_target(value: str) -> tuple[float, float]:
    fields = value.split(":")
    if len(fields) != 2:
        raise argparse.ArgumentTypeError("target must be PIM_REJECT:PI0_REJECT")
    target = tuple(float(field) for field in fields)
    if any(item <= 0.0 or item >= 1.0 for item in target):
        raise argparse.ArgumentTypeError("targets must be between zero and one")
    return target  # type: ignore[return-value]


def synthetic_self_test() -> None:
    rng = np.random.default_rng(6302026)
    parts = {
        "electron": rng.normal(0.0, 1.0, size=(400, len(BGO_FEATURES))),
        "pim": rng.normal(1.0, 1.0, size=(400, len(BGO_FEATURES))),
        "pi0": rng.normal(-0.8, 1.4, size=(400, len(BGO_FEATURES))),
    }
    for method in METHODS:
        model = fit_method(method, parts)
        scores = {
            species: score_method(method, model, values)
            for species, values in parts.items()
        }
        if method == "pairwise_qda":
            threshold = choose_pairwise_thresholds(scores, 0.80, 0.80)
            result = metrics_pairwise(scores, threshold)
        else:
            threshold = choose_threshold(scores, 0.80, 0.80)
            result = metrics(scores, threshold)
        if result["pim_reject"] < 0.80 or result["pi0_reject"] < 0.80:
            raise AssertionError((method, result))
        if not all(np.all(np.isfinite(value)) for value in scores.values()):
            raise AssertionError(f"non-finite score for {method}")
    print("bgo_compare_v1 synthetic self-test: OK")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ball-electron", type=Path)
    parser.add_argument("--ball-pim", type=Path)
    parser.add_argument("--ball-pi0", type=Path)
    parser.add_argument("--egg-electron", type=Path)
    parser.add_argument("--egg-pim", type=Path)
    parser.add_argument("--egg-pi0", type=Path)
    parser.add_argument("--ball-label", default="current_ball")
    parser.add_argument("--egg-label", default="bgoegg_envelope")
    parser.add_argument(
        "--feature-group", choices=sorted(FEATURE_GROUPS), default="bgo"
    )
    parser.add_argument(
        "--target",
        action="append",
        type=parse_target,
        default=None,
        help="validation targets PIM_REJECT:PI0_REJECT; repeatable",
    )
    parser.add_argument("--output", type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        synthetic_self_test()
        return
    required = (
        "ball_electron", "ball_pim", "ball_pi0", "egg_electron", "egg_pim",
        "egg_pi0", "output",
    )
    missing = [name for name in required if getattr(args, name) is None]
    if missing:
        parser.error("missing required arguments: " + ", ".join(missing))
    targets = args.target or [(0.93, 0.98), (0.94, 0.98), (0.94, 0.99)]
    feature_names = FEATURE_GROUPS[args.feature_group]
    methods = METHODS if args.feature_group != "all" else METHODS[:-1]
    ball = load_geometry(
        args.ball_electron,
        args.ball_pim,
        args.ball_pi0,
        args.ball_label,
        feature_names,
    )
    egg = load_geometry(
        args.egg_electron,
        args.egg_pim,
        args.egg_pi0,
        args.egg_label,
        feature_names,
    )
    result = {
        "status": "complete",
        "analysis": "bgo_compare_v1",
        "feature_group": args.feature_group,
        "classifier_inputs": (
            feature_names if feature_names is not None
            else "all non-eventID columns in each geometry"
        ),
        "truth_usage": "sample label only; no pdg/pid/truth path/true vertex",
        "split": {
            "fit": "eventID mod 4 == 0",
            "validation": "eventID mod 4 == 2; threshold and method selection",
            "test": "odd eventID; final measurement only",
        },
        "targets": [
            {"pim_reject": target[0], "pi0_reject": target[1]}
            for target in targets
        ],
        "methods": list(methods),
        "geometries": {
            "ball": analyze_geometry(ball, targets, methods),
            "egg": analyze_geometry(egg, targets, methods),
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    summary = {}
    for geometry_name, geometry in result["geometries"].items():
        summary[geometry_name] = {
            key: {
                "selected_method": value["selected_method"],
                **value["selected_test"],
            }
            for key, value in geometry["targets"].items()
        }
    print(json.dumps(summary, indent=2))
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
