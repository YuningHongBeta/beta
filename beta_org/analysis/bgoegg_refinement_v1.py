#!/usr/bin/env python3
"""Analyze the combined BGOegg end-ring and z-offset refinement scan."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from bgo_compare_v1 import FEATURE_GROUPS, METHODS, analyze_geometry
from bgoegg_study_v1 import (
    DEFAULT_CACHE,
    DEFAULT_OUTPUT as DEFAULT_BASE_RESULT,
    PROJECT,
    extract_geometry,
    load_triplet,
    load_valid_state,
    roots_from_state,
    source_manifest,
)


DEFAULT_STATE = Path(
    "/gpfs/home/had/yhong/k18-analyzer-e63/tmp/"
    "beta_production_t1_worktree/beta_org/runmanager/.state/"
    "bgoegg_refinement_t1.json"
)
DEFAULT_OUTPUT = (
    PROJECT / "analysis" / "results" / "bgoegg_v1" /
    "bgoegg_refinement_v1.json"
)


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
    parser.add_argument("--state", type=Path, default=DEFAULT_STATE)
    parser.add_argument("--base-result", type=Path, default=DEFAULT_BASE_RESULT)
    parser.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    for executable in (args.feature_extractor, args.pattern_extractor):
        if not executable.is_file():
            parser.error(f"missing extractor: {executable}")
    base = json.loads(args.base_result.read_text(encoding="utf-8"))
    if base.get("status") != "complete":
        parser.error(f"incomplete base result: {args.base_result}")
    definition = base["comparison_definition"]
    target = (
        float(definition["target_pim_reject"]),
        float(definition["target_pi0_reject"]),
    )

    state = load_valid_state(args.state)
    roots = roots_from_state(state, "refinement")
    standard = {}
    pattern = {}
    for label, triplet in roots.items():
        standard_files = extract_geometry(
            args.feature_extractor, triplet, args.cache, label, "standard", 1.0
        )
        pattern_files = extract_geometry(
            args.pattern_extractor, triplet, args.cache, label, "pattern", 1.0
        )
        standard[label] = analyze_geometry(
            load_triplet(standard_files, label, FEATURE_GROUPS["bgo"]),
            [target],
            METHODS,
        )
        pattern[label] = analyze_geometry(
            load_triplet(pattern_files, label, None),
            [target],
            METHODS[:-1],
        )

    result = {
        "status": "complete",
        "analysis": "bgoegg_refinement_v1",
        "base_result": str(args.base_result.resolve()),
        "comparison_definition": definition,
        "truth_usage": "species sample label only; no event truth classifier input",
        "state": {
            "path": str(args.state.resolve()),
            "schema": state["schema"],
            "provenance": state["provenance"],
        },
        "source_roots": source_manifest(roots),
        "standard_bgo_hardware": standard,
        "pattern_bgo_hardware": pattern,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(args.output)


if __name__ == "__main__":
    main()
