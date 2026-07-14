#!/usr/bin/env python3
"""Separate PC gamma interception and intrinsic/effective veto response."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import yaml

from bgo_features_v2 import load_bgo2


PROJECT = Path(__file__).resolve().parent.parent
TAG = "bgoegg_pc_aperture_s6402026"
DEFAULT_ROOT = PROJECT / "build_bgoegg_frustum" / "output" / "scan"
DEFAULT_CACHE = PROJECT / "tmp" / TAG
DEFAULT_MANIFEST = PROJECT / "runmanager" / "examples" / f"{TAG}.yml"
DEFAULT_OUTPUT = (
    PROJECT / "analysis" / "results" / "bgoegg_v1" / f"{TAG}.json"
)
DEFAULT_EXTRACTOR = PROJECT / "analysis" / "bgo_extract_features_v2"
PC_THRESHOLDS_MEV = (0.0, 0.1, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0)


def wilson(p: float, n: int):
    if n <= 0:
        return None
    z = 1.959963984540054
    denominator = 1.0 + z * z / n
    center = (p + z * z / (2.0 * n)) / denominator
    half = z * math.sqrt(
        p * (1.0 - p) / n + z * z / (4.0 * n * n)
    ) / denominator
    return [float(center - half), float(center + half)]


def measured(mask: np.ndarray, condition=None) -> dict:
    if condition is None:
        condition = np.ones(len(mask), dtype=bool)
    denominator = int(np.count_nonzero(condition))
    numerator = int(np.count_nonzero(mask & condition))
    value = None if denominator == 0 else numerator / denominator
    return {
        "value": value, "numerator": numerator, "denominator": denominator,
        "wilson_95ci": None if value is None else wilson(value, denominator),
    }


def ensure_features(root: Path, output: Path, extractor: Path) -> None:
    if output.is_file() and Path(str(output) + ".json").is_file():
        return
    output.parent.mkdir(parents=True, exist_ok=True)
    import subprocess
    subprocess.run([str(extractor), str(root), str(output), "3"], check=True)


def masks(data: np.ndarray, names: list[str]) -> dict[str, np.ndarray]:
    col = {name: index for index, name in enumerate(names)}
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
    return {"all": np.ones(len(data), dtype=bool), "bgo": bgo,
            "pre_pc": bgo & th & tlc & dz}


def hit(data: np.ndarray, index: int, threshold: float) -> np.ndarray:
    return data[:, index] > 0.0 if threshold == 0.0 else data[:, index] >= threshold


def evaluate_geometry(tag: str, name: str, design: dict, root_dir: Path,
                      cache: Path, extractor: Path) -> dict:
    samples, names = {}, None
    roots = {}
    for species in ("e", "pi0"):
        root = root_dir / f"{tag}_{name}_{species}.root"
        feature = cache / name / f"{species}.bgo2"
        ensure_features(root, feature, extractor)
        data, meta, manifest = load_bgo2(feature, copy=False)
        if not math.isclose(float(meta["threshold_MeV"]), 3.0):
            raise ValueError("BGO threshold is not 3 MeV")
        if names is None:
            names = list(manifest["features"])
        elif names != list(manifest["features"]):
            raise ValueError("feature order mismatch")
        samples[species] = data
        roots[species] = manifest["input"]
    assert names is not None
    col = {feature: index for index, feature in enumerate(names)}
    required = {
        "pcGammaN", "pcGammaDownN", "pcGammaUpN", "pcSumE_MeV",
        "pcDownE_MeV", "pcUpE_MeV",
    }
    if not required.issubset(col):
        raise ValueError(f"gamma-entrance features absent for {name}")
    selection = {species: masks(data, names) for species, data in samples.items()}
    intercept = {
        species: data[:, col["pcGammaN"]] > 0.0
        for species, data in samples.items()
    }
    intercept_up = {
        species: data[:, col["pcGammaUpN"]] > 0.0
        for species, data in samples.items()
    }
    intercept_down = {
        species: data[:, col["pcGammaDownN"]] > 0.0
        for species, data in samples.items()
    }
    rows = []
    for threshold in PC_THRESHOLDS_MEV:
        hits = {
            species: hit(data, col["pcSumE_MeV"], threshold)
            for species, data in samples.items()
        }
        row = {
            "pc_plastic_threshold_MeV": threshold,
            "hit_rule": "EdepPC>0" if threshold == 0.0 else (
                f"EdepPC >= {threshold:g} MeV"),
            "all_generated": {}, "given_bgo_ncluster1": {},
            "given_full_pre_pc": {},
        }
        for species in ("e", "pi0"):
            for label, condition_name in (
                    ("all_generated", "all"),
                    ("given_bgo_ncluster1", "bgo"),
                    ("given_full_pre_pc", "pre_pc")):
                condition = selection[species][condition_name]
                row[label][species] = {
                    "gamma_intercept": measured(intercept[species], condition),
                    "effective_pc_hit": measured(hits[species], condition),
                    "intrinsic_hit_given_gamma_intercept": measured(
                        hits[species], condition & intercept[species]),
                    "pc_hit_without_recorded_gamma_intercept": measured(
                        hits[species] & ~intercept[species], condition),
                }
        rows.append(row)
    return {
        "design": design,
        "source_root_files": roots,
        "selection_counts": {
            species: {
                key: int(np.count_nonzero(value))
                for key, value in selection[species].items()
            } for species in ("e", "pi0")
        },
        "gamma_intercept_side_all_generated": {
            species: {
                "any": measured(intercept[species]),
                "upstream": measured(intercept_up[species]),
                "downstream": measured(intercept_down[species]),
            } for species in ("e", "pi0")
        },
        "operating_points": rows,
    }


def point_at(result: dict, threshold: float) -> dict:
    return next(row for row in result["operating_points"]
                if row["pc_plastic_threshold_MeV"] == threshold)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root-dir", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--extractor", type=Path, default=DEFAULT_EXTRACTOR)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    manifest = yaml.safe_load(args.manifest.read_text(encoding="utf-8"))
    tag = str(manifest["tag"])
    designs = {item["name"]: item for item in manifest["geometries"]}
    geometries = {
        name: evaluate_geometry(tag, name, design, args.root_dir.resolve(),
                                args.cache.resolve(), args.extractor.resolve())
        for name, design in designs.items()
    }
    threshold = 0.5
    screen = []
    for name, geometry in geometries.items():
        point = point_at(geometry, threshold)
        screen.append({
            "geometry": name,
            "pc_plastic_threshold_MeV": threshold,
            "pi0_given_bgo": point["given_bgo_ncluster1"]["pi0"],
            "pi0_given_full_pre_pc": point["given_full_pre_pc"]["pi0"],
            "electron_given_bgo": point["given_bgo_ncluster1"]["e"],
        })
    selected = max(
        screen,
        key=lambda item: (
            item["pi0_given_bgo"]["effective_pc_hit"]["value"] or -1.0,
            -(item["electron_given_bgo"]["effective_pc_hit"]["value"] or 0.0),
        ),
    )
    result = {
        "status": "complete",
        "analysis": "bgoegg_pc_aperture_evaluate",
        "tag": tag,
        "bgo_threshold_MeV": 3.0,
        "pc_policy": "hard stand-alone PC hit veto; no joint classifier",
        "screen_threshold_MeV": threshold,
        "selected_ideal_upper_bound": selected,
        "screening_table": screen,
        "geometry_warning": (
            "full-aperture axisymmetric candidates cover charged-particle "
            "corridors and are ideal upper bounds, not installable designs"
        ),
        "thesis_benchmark": {
            "ideal_pc_veto_given_remaining_pi0": 0.70,
            "kamada_realistic_100MeV_gamma_detection_max": 0.64,
        },
        "definitions": {
            "gamma_intercept": (
                "at least one unique gamma track enters the first PC Pb layer"
            ),
            "intrinsic": "PC plastic hit conditional on gamma_intercept",
            "effective": "PC plastic hit including geometry and conversion/response",
        },
        "geometries": geometries,
        "manifest": str(args.manifest.resolve()),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"selected_ideal_upper_bound": selected}, indent=2))
    print(args.output)


if __name__ == "__main__":
    main()
