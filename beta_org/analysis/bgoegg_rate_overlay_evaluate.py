#!/usr/bin/env python3
"""Apply the frozen exact-BGOegg T3 classifier to K1.8 rate overlays.

The clean seed-6302026 sample alone defines the model and score rectangle.
No beam-contaminated sample is used for fitting, feature selection, or cut
selection.  The seed-8502026 clean/overlay samples are evaluation-only.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
from pathlib import Path

import numpy as np

from bgo_compare_v1 import choose_pairwise_thresholds, metrics_pairwise, split_masks
from bgo_threshold3_analysis_v1 import (
    HURDLE_FEATURE_SETS,
    fit_hurdle,
    hurdle_feature_values,
    hurdle_scores,
)
from bgoegg_frustum31_threshold3_scan import ROBUST_TARGET


PROJECT = Path(__file__).resolve().parent.parent
DEFAULT_BUILD = PROJECT / "build_rateoverlay" / "output" / "scan"
DEFAULT_TRAINING = (
    PROJECT / "tmp" / "bgoegg_frustum31_zscan_s6302026" / "extended31_zm10"
)
DEFAULT_CACHE = PROJECT / "tmp" / "bgoegg_rate_overlay_s8502026"
DEFAULT_OUTPUT = (
    PROJECT / "analysis" / "results" / "bgoegg_v1" /
    "bgoegg_rate_overlay_s8502026.json"
)
EXTRACTOR = PROJECT / "analysis" / "bgo_extract_features_v2"
FEATURE_SET = "energy_rms_isolation"
SPECIES = ("electron", "pim", "pi0")
PRIMARY = {"electron": "e", "pim": "pim", "pi0": "pi0"}
CAMPAIGNS = {
    "0.30": "hong_k18_bgoegg31_rate_low_lam0p30_s8502026",
    "0.595": "hong_k18_bgoegg31_rate_nominal_s8502026",
    "1.25": "hong_k18_bgoegg31_rate_high_lam1p25_s8502026",
}
GEOMETRY = "extended31_zm10"
GATE_TAG = "hg10p0ns"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def root_path(build: Path, tag: str, primary: str) -> Path:
    return build / f"{tag}_{GEOMETRY}_{GATE_TAG}_{primary}.root"


def extract_one(source: Path, output: Path, force: bool) -> None:
    sidecar = Path(str(output) + ".json")
    if output.is_file() and sidecar.is_file() and not force:
        return
    if not source.is_file():
        raise FileNotFoundError(source)
    output.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run([str(EXTRACTOR), str(source), str(output), "3"], check=True)


def prepare_inputs(build: Path, cache: Path, force: bool) -> tuple[Path, dict[str, Path]]:
    clean = cache / "clean"
    nominal_tag = CAMPAIGNS["0.595"]
    for species in SPECIES:
        primary = PRIMARY[species]
        extract_one(root_path(build, nominal_tag, primary), clean / f"{primary}.bgo2", force)
    overlays = {}
    for rate, tag in CAMPAIGNS.items():
        directory = cache / f"rate_{rate.replace('.', 'p')}"
        overlays[rate] = directory
        for species in SPECIES:
            primary = PRIMARY[species]
            extract_one(
                root_path(build, tag, f"{primary}_beam"),
                directory / f"{primary}.bgo2",
                force,
            )
    extract_one(root_path(build, nominal_tag, "beam"), cache / "beam.bgo2", force)
    return clean, overlays


def score_geometry(directory: Path, label: str, model: dict) -> tuple[dict, dict]:
    geometry = hurdle_feature_values(
        directory, label, HURDLE_FEATURE_SETS[FEATURE_SET]
    )
    scores = {
        species: hurdle_scores(
            model, sample, np.ones(len(sample["event_id"]), dtype=bool)
        )
        for species, sample in geometry["samples"].items()
    }
    return geometry, scores


def accepted(scores: np.ndarray, cuts: tuple[float, float]) -> np.ndarray:
    return (scores[:, 0] < cuts[0]) & (scores[:, 1] < cuts[1])


def paired_change(clean_geometry: dict, clean_scores: dict, overlay_geometry: dict,
                  overlay_scores: dict, cuts: tuple[float, float]) -> dict:
    result = {}
    for species in SPECIES:
        clean_sample = clean_geometry["samples"][species]
        overlay_sample = overlay_geometry["samples"][species]
        clean_order = np.argsort(clean_sample["event_id"])
        overlay_order = np.argsort(overlay_sample["event_id"])
        clean_ids = clean_sample["event_id"][clean_order]
        overlay_ids = overlay_sample["event_id"][overlay_order]
        if not np.array_equal(clean_ids, overlay_ids):
            raise ValueError(f"event IDs do not pair for {species}")
        before = accepted(clean_scores[species], cuts)[clean_order]
        after = accepted(overlay_scores[species], cuts)[overlay_order]
        pass_to_fail = int(np.count_nonzero(before & ~after))
        fail_to_pass = int(np.count_nonzero(~before & after))
        nvalue = len(before)
        delta = (fail_to_pass - pass_to_fail) / nvalue
        discordant = pass_to_fail + fail_to_pass
        variance = max(
            0.0,
            (discordant - (fail_to_pass - pass_to_fail) ** 2 / nvalue)
            / (nvalue * nvalue),
        )
        half = 1.959963984540054 * math.sqrt(variance)
        result[species] = {
            "total": nvalue,
            "pass_to_fail": pass_to_fail,
            "fail_to_pass": fail_to_pass,
            "pass_both": int(np.count_nonzero(before & after)),
            "fail_both": int(np.count_nonzero(~before & ~after)),
            "overlay_minus_clean": delta,
            "paired_normal_95ci": [delta - half, delta + half],
        }
    return result


def discrete_distribution(geometry: dict, feature: str) -> dict:
    result = {}
    for species, sample in geometry["samples"].items():
        values = np.asarray(sample["raw"][:, sample["names"].index(feature)], dtype=int)
        counts = [int(np.count_nonzero(values == value)) for value in range(5)]
        counts.append(int(np.count_nonzero(values >= 5)))
        result[species] = {"bins": ["0", "1", "2", "3", "4", "5+"], "counts": counts}
    return result


def score_histograms(all_scores: dict[str, dict[str, np.ndarray]]) -> dict:
    result = {"definition": "under/overflow are included in the edge bins"}
    for axis, name in enumerate(("pim_minus_e_log_likelihood", "pi0_minus_e_log_likelihood")):
        pooled = np.concatenate([
            values[:, axis]
            for condition in all_scores.values()
            for values in condition.values()
        ])
        low, high = np.quantile(pooled, [0.001, 0.999])
        if not np.isfinite(low) or not np.isfinite(high) or low >= high:
            low, high = float(np.min(pooled)), float(np.max(pooled) + 1.0)
        edges = np.linspace(low, high, 81)
        axis_result = {"edges": edges.tolist(), "conditions": {}}
        for condition, species_scores in all_scores.items():
            axis_result["conditions"][condition] = {}
            for species, values in species_scores.items():
                clipped = np.clip(values[:, axis], edges[0], edges[-1])
                axis_result["conditions"][condition][species] = np.histogram(
                    clipped, bins=edges
                )[0].astype(int).tolist()
        result[name] = axis_result
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--build", type=Path, default=DEFAULT_BUILD)
    parser.add_argument("--training", type=Path, default=DEFAULT_TRAINING)
    parser.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--force-extract", action="store_true")
    args = parser.parse_args()

    clean_dir, overlay_dirs = prepare_inputs(
        args.build.resolve(), args.cache.resolve(), args.force_extract
    )
    training = hurdle_feature_values(
        args.training.resolve(), "seed6302026-training",
        HURDLE_FEATURE_SETS[FEATURE_SET],
    )
    model = fit_hurdle(training)
    validation_scores = {
        species: hurdle_scores(
            model, sample, split_masks(sample["event_id"])["validation"]
        )
        for species, sample in training["samples"].items()
    }
    cuts = choose_pairwise_thresholds(validation_scores, *ROBUST_TARGET)

    clean_geometry, clean_scores = score_geometry(clean_dir, "clean", model)
    overlay_geometry = {}
    overlay_scores = {}
    for rate, directory in overlay_dirs.items():
        overlay_geometry[rate], overlay_scores[rate] = score_geometry(
            directory, f"overlay-rate-{rate}", model
        )

    conditions = {"clean": clean_scores}
    conditions.update({f"rate_{rate}": value for rate, value in overlay_scores.items()})
    inputs = sorted(args.build.resolve().glob(
        "hong_k18_bgoegg31_rate_*_s8502026_extended31_zm10_hg10p0ns_*.root"
    ))
    output = {
        "status": "complete",
        "analysis": "bgoegg_rate_overlay_evaluate",
        "definition": (
            "frozen seed-6302026 exact-BGOegg hurdle model and robust validation "
            "cuts applied without refitting to seed-8502026 clean and K1.8 Poisson overlays"
        ),
        "geometry": {
            "name": GEOMETRY, "rings": 31, "sectors": 60,
            "center_z_cm": -10, "per_cell_threshold_MeV": 3.0,
        },
        "beam_model": {
            "particle": "K-", "momentum_GeV_c": 1.5,
            "poisson_means_per_1us": [0.30, 0.595, 1.25],
            "profile": "independent truncated Gaussian",
            "sigma_x_mm": 24, "sigma_y_mm": 5,
            "max_abs_x_mm": 60, "max_abs_y_mm": 20,
        },
        "analysis_contract": {
            "feature_set": FEATURE_SET,
            "features": HURDLE_FEATURE_SETS[FEATURE_SET],
            "model": "nHit-binned diagonal Gaussian likelihood",
            "fit": "seed6302026 eventID mod 4 == 0",
            "cut_selection": "seed6302026 eventID mod 4 == 2",
            "target": {"pim_reject": ROBUST_TARGET[0], "pi0_reject": ROBUST_TARGET[1]},
            "thresholds": [float(value) for value in cuts],
            "training_validation": metrics_pairwise(validation_scores, cuts),
            "contaminated_sample_used_for_training_or_cut_selection": False,
        },
        "metrics": {
            "clean": metrics_pairwise(clean_scores, cuts),
            **{
                f"rate_{rate}": metrics_pairwise(scores, cuts)
                for rate, scores in overlay_scores.items()
            },
        },
        "paired_changes": {
            f"rate_{rate}": paired_change(
                clean_geometry, clean_scores, overlay_geometry[rate], scores, cuts
            )
            for rate, scores in overlay_scores.items()
        },
        "ncluster4_distributions": {
            "clean": discrete_distribution(clean_geometry, "nCl4"),
            **{
                f"rate_{rate}": discrete_distribution(geometry, "nCl4")
                for rate, geometry in overlay_geometry.items()
            },
        },
        "nhit_distributions": {
            "clean": discrete_distribution(clean_geometry, "nHit"),
            **{
                f"rate_{rate}": discrete_distribution(geometry, "nHit")
                for rate, geometry in overlay_geometry.items()
            },
        },
        "score_histograms": score_histograms(conditions),
        "inputs": [
            {"path": str(path), "sha256": sha256(path)} for path in inputs
        ],
        "limitations": [
            "K1.8 rate and target-plane profile are proxies, not a digitized measured spill model",
            "BGO response is a 1 us rectangular Geant4-step gate without waveform shaping, dead time, or electronics noise",
            "13C target transverse size/material density and free-Lambda signal are transport proxies",
            "support, PMT, gain spread, energy resolution, dead channels, and S-2S material are absent",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output["metrics"], indent=2))
    print(args.output)


if __name__ == "__main__":
    main()
