#!/usr/bin/env python3
"""Select a 3 MeV analysis and z position for the exact 31-ring BGOegg.

The geometry/method choice uses only eventID mod 4 == 2.  Odd event IDs are
reported as a candidate evaluation and must not be used to revise the choice.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from bgo_compare_v1 import FEATURE_GROUPS, METHODS, analyze_geometry
from bgo_threshold3_analysis_v1 import (
    TARGET,
    compact_scan,
    hurdle_scan,
    load_triplet,
    ncluster_distribution,
    selected_summary,
)


PROJECT = Path(__file__).resolve().parent.parent
DEFAULT_INPUT = PROJECT / "tmp" / "bgoegg_frustum31_zscan_s6302026"
DEFAULT_OUTPUT = (
    PROJECT
    / "analysis"
    / "results"
    / "bgoegg_v1"
    / "bgoegg_frustum31_threshold3_scan_s6302026.json"
)
GEOMETRIES = (
    "nominal22_z0",
    "extended31_zm10",
    "extended31_zm5",
    "extended31_z0",
    "extended31_zp5",
    "extended31_zp10",
)
REFERENCE = {
    "label": "current ball at 1 MeV, seed 6302026",
    "electron_keep": 0.60156,
    "pim_reject": TARGET[0],
    "pi0_reject": TARGET[1],
}


def candidate_rows(label: str, result: dict) -> list[dict]:
    standard = result["standard"]
    compact = result["compact"]
    hurdle = result["hurdle"]
    return [
        {
            "geometry": label,
            "family": "standard_full",
            "feature_set": "bgo",
            "method": standard["selected_method"],
            "validation": standard["validation"],
            "candidate_evaluation": standard["candidate_evaluation"],
        },
        {
            "geometry": label,
            "family": "compact",
            "feature_set": compact["selected_feature_set"],
            "method": compact["selected"]["selected_method"],
            "validation": compact["selected"]["validation"],
            "candidate_evaluation": compact["selected"]["candidate_evaluation"],
        },
        {
            "geometry": label,
            "family": "hurdle",
            "feature_set": hurdle["selected_feature_set"],
            "method": "nHit-binned diagonal Gaussian likelihood",
            "validation": hurdle["selected"]["validation"],
            "candidate_evaluation": hurdle["selected"]["candidate_evaluation"],
        },
        {
            "geometry": label,
            "family": "hurdle_guarded",
            "feature_set": hurdle["selected_feature_set"],
            "method": "nHit-binned diagonal Gaussian likelihood",
            "validation": hurdle["guarded_operating_point"]["validation"],
            "candidate_evaluation": hurdle["guarded_operating_point"][
                "candidate_evaluation"
            ],
        },
    ]


def analyze(directory: Path, label: str) -> dict:
    standard_result = analyze_geometry(
        load_triplet(directory, label + ":standard", FEATURE_GROUPS["bgo"]),
        [TARGET],
        METHODS,
    )
    return {
        "standard": selected_summary(standard_result),
        "compact": compact_scan(directory, label),
        "hurdle": hurdle_scan(directory, label),
        "ncluster_full_sample": ncluster_distribution(
            directory, label + ":ncluster"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    args.input = args.input.resolve()

    geometries = {
        label: analyze(args.input / label, label) for label in GEOMETRIES
    }
    rows = [
        row
        for label, result in geometries.items()
        for row in candidate_rows(label, result)
    ]
    # This is the only geometry/method selection.  Odd-event values below are
    # retained solely to estimate the already selected candidate's performance.
    selected = max(rows, key=lambda row: row["validation"]["e_keep"])
    evaluation = selected["candidate_evaluation"]
    result = {
        "status": "complete",
        "analysis": "bgoegg_frustum31_threshold3_scan",
        "input": str(args.input),
        "threshold_MeV": 3.0,
        "threshold_policy": (
            "hard per-segment cut; deposits below 3 MeV are not classifier inputs"
        ),
        "reference": REFERENCE,
        "split": {
            "fit": "eventID mod 4 == 0",
            "validation": "eventID mod 4 == 2; geometry/method/feature/cut selection",
            "candidate_evaluation": (
                "odd eventID; inspected only after the single selection"
            ),
        },
        "selection_rule": (
            "maximum validation electron retention among all declared geometry, "
            "method-family, and feature-set candidates"
        ),
        "selected": selected,
        "selected_gap_to_ball_T1": {
            "electron_keep_percentage_point": 100.0
            * (evaluation["e_keep"] - REFERENCE["electron_keep"]),
            "pim_reject_percentage_point": 100.0
            * (evaluation["pim_reject"] - REFERENCE["pim_reject"]),
            "pi0_reject_percentage_point": 100.0
            * (evaluation["pi0_reject"] - REFERENCE["pi0_reject"]),
        },
        "candidate_table": rows,
        "geometries": geometries,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "selected": selected,
        "selected_gap_to_ball_T1": result["selected_gap_to_ball_T1"],
    }, indent=2))
    print(args.output)


if __name__ == "__main__":
    main()
