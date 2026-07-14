#!/usr/bin/env python3
"""Run the reproducible BGO-only BGOegg analysis and hardware comparison."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess

from bgo_compare_v1 import (
    FEATURE_GROUPS,
    METHODS,
    analyze_geometry,
    cluster_baseline,
    load_geometry,
)
import numpy as np


PROJECT = Path(__file__).resolve().parent.parent
REPOSITORY = PROJECT.parent
DEFAULT_BALL_DIR = PROJECT / "build_thtiming" / "output" / "scan"
DEFAULT_PRODUCTION_STATE = Path(
    "/gpfs/home/had/yhong/k18-analyzer-e63/tmp/beta_production_t1_worktree/"
    "beta_org/runmanager/.state/bgoegg_prod_t1.json"
)
DEFAULT_CACHE = PROJECT / "tmp" / "bgoegg_study_v1"
DEFAULT_OUTPUT = (
    PROJECT / "analysis" / "results" / "bgoegg_v1" /
    "bgoegg_study_v1.json"
)
THRESHOLDS = (0.25, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0)
SPECIES = ("e", "pim", "pi0")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def threshold_tag(threshold: float) -> str:
    return f"T{threshold:g}".replace(".", "p")


def load_valid_state(path: Path) -> dict:
    state = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(state.get("jobs"), list) or not state["jobs"]:
        raise RuntimeError(f"invalid runmanager state: {path}")
    invalid = [
        job["key"] for job in state["jobs"]
        if job.get("validation", {}).get("status") != "valid"
    ]
    if invalid:
        raise RuntimeError(f"unvalidated jobs in {path}: {invalid}")
    for job in state["jobs"]:
        if not Path(job["output_root"]).is_file():
            raise RuntimeError(f"missing validated ROOT: {job['output_root']}")
    return state


def roots_from_state(state: dict, state_label: str) -> dict[str, dict[str, Path]]:
    result: dict[str, dict[str, Path]] = {}
    for job in state["jobs"]:
        geometry = job["geometry"]["name"]
        result.setdefault(f"{state_label}:{geometry}", {})[job["primary"]] = Path(
            job["output_root"]
        )
    for label, roots in result.items():
        if set(roots) != set(SPECIES):
            raise RuntimeError(f"{label}: expected {SPECIES}, got {sorted(roots)}")
    return result


def extract_one(
    executable: Path,
    source: Path,
    destination: Path,
    threshold: float,
) -> None:
    sidecar = Path(str(destination) + ".json")
    if destination.is_file() and sidecar.is_file():
        manifest = json.loads(sidecar.read_text(encoding="utf-8"))
        if (
            Path(manifest.get("input", "")).resolve() == source.resolve()
            and float(manifest.get("threshold_MeV", -1.0)) == threshold
            and int(manifest.get("nrow", -1)) == 100_000
        ):
            return
    destination.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [str(executable), str(source), str(destination), str(threshold)],
        check=True,
    )


def extract_geometry(
    executable: Path,
    roots: dict[str, Path],
    cache: Path,
    label: str,
    kind: str,
    threshold: float,
) -> dict[str, Path]:
    output = {}
    safe_label = label.replace(":", "__")
    for species in SPECIES:
        path = cache / kind / safe_label / (
            f"{species}_{threshold_tag(threshold)}.bgo2"
        )
        extract_one(executable, roots[species], path, threshold)
        output[species] = path
    return output


def load_triplet(
    files: dict[str, Path], label: str, feature_names: list[str] | None
) -> dict:
    return load_geometry(
        files["e"], files["pim"], files["pi0"], label, feature_names
    )


def combine_geometries(label: str, standard: dict, pattern: dict) -> dict:
    """Join summary and cell-pattern inputs by eventID for one ablation."""
    samples = {}
    for species in ("electron", "pim", "pi0"):
        first = standard["samples"][species]
        second = pattern["samples"][species]
        if not np.array_equal(first["event_id"], second["event_id"]):
            raise RuntimeError(f"{label}: standard/pattern eventID mismatch for {species}")
        combined = dict(first)
        combined["x"] = np.column_stack([first["x"], second["x"]])
        combined["classifier_inputs"] = (
            first["classifier_inputs"] + second["classifier_inputs"]
        )
        samples[species] = combined
    return {
        "label": label,
        "samples": samples,
        "metadata": {"combined_standard_and_pattern": True},
        "feature_names": standard["feature_names"] + pattern["feature_names"],
    }


def source_manifest(roots: dict[str, dict[str, Path]]) -> dict:
    return {
        label: {
            species: {"path": str(path), "sha256": sha256(path)}
            for species, path in triplet.items()
        }
        for label, triplet in roots.items()
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--feature-extractor", type=Path,
        default=PROJECT / "analysis" / "bgo_extract_features_v2",
    )
    parser.add_argument(
        "--pattern-extractor", type=Path,
        default=PROJECT / "analysis" / "bgo_extract_pattern_v1",
    )
    parser.add_argument("--ball-dir", type=Path, default=DEFAULT_BALL_DIR)
    parser.add_argument(
        "--production-state", type=Path, default=DEFAULT_PRODUCTION_STATE
    )
    parser.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    for executable in (args.feature_extractor, args.pattern_extractor):
        if not executable.is_file():
            parser.error(f"missing extractor: {executable}")

    states = {"production": load_valid_state(args.production_state)}
    ball_roots = {
        species: args.ball_dir / f"timing_a20x40_{species}.root"
        for species in SPECIES
    }
    if any(not path.is_file() for path in ball_roots.values()):
        parser.error(f"missing authoritative ball ROOT under {args.ball_dir}")
    roots = {"ball:a20x40": ball_roots}
    for state_label, state in states.items():
        roots.update(roots_from_state(state, state_label))

    standard_files: dict[str, dict[str, Path]] = {}
    pattern_files: dict[str, dict[str, Path]] = {}
    for label, triplet in roots.items():
        standard_files[label] = extract_geometry(
            args.feature_extractor, triplet, args.cache, label, "standard", 1.0
        )
        pattern_files[label] = extract_geometry(
            args.pattern_extractor, triplet, args.cache, label, "pattern", 1.0
        )

    ball_standard = load_triplet(
        standard_files["ball:a20x40"], "ball:a20x40", FEATURE_GROUPS["bgo"]
    )
    ball_cluster = cluster_baseline(ball_standard)
    if ball_cluster is None:
        raise RuntimeError("ball nCl4 baseline is unavailable")
    target = (
        ball_cluster["test"]["pim_reject"],
        ball_cluster["test"]["pi0_reject"],
    )
    targets = [target]

    standard_hardware = {
        "ball:a20x40": analyze_geometry(ball_standard, targets, METHODS)
    }
    for label, files in standard_files.items():
        if label == "ball:a20x40":
            continue
        geometry = load_triplet(files, label, FEATURE_GROUPS["bgo"])
        standard_hardware[label] = analyze_geometry(geometry, targets, METHODS)

    ball_pattern = load_triplet(
        pattern_files["ball:a20x40"], "ball:a20x40", None
    )
    pattern_methods = METHODS[:-1]
    pattern_hardware = {
        "ball:a20x40": analyze_geometry(ball_pattern, targets, pattern_methods)
    }
    for label, files in pattern_files.items():
        if label == "ball:a20x40":
            continue
        geometry = load_triplet(files, label, None)
        pattern_hardware[label] = analyze_geometry(
            geometry, targets, pattern_methods
        )

    combined_baseline = {}
    for label in ("ball:a20x40", "production:egg_none"):
        standard_geometry = load_triplet(
            standard_files[label], label + ":standard", FEATURE_GROUPS["bgo"]
        )
        pattern_geometry = load_triplet(
            pattern_files[label], label + ":pattern", None
        )
        combined_baseline[label] = analyze_geometry(
            combine_geometries(
                label + ":combined", standard_geometry, pattern_geometry
            ),
            targets,
            pattern_methods,
        )

    threshold_scan = {}
    baseline_label = "production:egg_none"
    for threshold in THRESHOLDS:
        files_by_label = {}
        for label in ("ball:a20x40", baseline_label):
            files_by_label[label] = extract_geometry(
                args.feature_extractor,
                roots[label],
                args.cache,
                label,
                "standard",
                threshold,
            )
        threshold_scan[threshold_tag(threshold)] = {
            label: analyze_geometry(
                load_triplet(files, label, FEATURE_GROUPS["bgo"]),
                targets,
                METHODS,
            )
            for label, files in files_by_label.items()
        }

    result = {
        "status": "complete",
        "analysis": "bgoegg_study_v1",
        "comparison_definition": {
            "reference": "authoritative current-ball equal-solid-angle 20x40",
            "reference_selection": "nCl4 == 1 at 1 MeV per-cell threshold",
            "target_pim_reject": target[0],
            "target_pi0_reject": target[1],
            "threshold_selection": "validation only",
            "method_selection": "maximum validation electron retention",
            "final_measurement": "odd eventID test subset",
        },
        "truth_usage": "species sample label only; no event truth classifier input",
        "extractors": {
            "standard": {
                "path": str(args.feature_extractor.resolve()),
                "sha256": sha256(args.feature_extractor),
            },
            "pattern": {
                "path": str(args.pattern_extractor.resolve()),
                "sha256": sha256(args.pattern_extractor),
            },
        },
        "states": {
            name: {
                "path": str(path.resolve()),
                "schema": states[name]["schema"],
                "provenance": states[name]["provenance"],
            }
            for name, path in (("production", args.production_state),)
        },
        "source_roots": source_manifest(roots),
        "ball_cluster_reference": ball_cluster,
        "threshold_scan": threshold_scan,
        "standard_bgo_hardware": standard_hardware,
        "pattern_bgo_hardware": pattern_hardware,
        "combined_standard_pattern_baseline": combined_baseline,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(args.output)


if __name__ == "__main__":
    main()
