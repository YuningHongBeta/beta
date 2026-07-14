#!/usr/bin/env python3
"""Select a 3 MeV analysis and z position for the exact 31-ring BGOegg.

The geometry/model/feature choice uses only eventID mod 4 == 2.  Odd event IDs
diagnose whether that discovery cut needs a fixed safety margin before the
independent seed is read; they never revise the geometry, model, or features.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from bgo_compare_v1 import FEATURE_GROUPS, METHODS, analyze_geometry
from bgo_threshold3_analysis_v1 import (
    HURDLE_FEATURE_SETS,
    TARGET,
    analyze_hurdle,
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
# Fixed after the exact-target discovery candidate missed the seed-6302026 odd
# background rates, and before any seed-7302026 confirmation output was read.
ROBUST_TARGET = (0.927, 0.9901)


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
        {
            "geometry": label,
            "family": "hurdle_robust",
            "feature_set": hurdle["selected_feature_set"],
            "method": "nHit-binned diagonal Gaussian likelihood",
            "validation": result["hurdle_robust"]["validation"],
            "candidate_evaluation": result["hurdle_robust"][
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
    hurdle = hurdle_scan(directory, label)
    hurdle_feature_set = hurdle["selected_feature_set"]
    return {
        "standard": selected_summary(standard_result),
        "compact": compact_scan(directory, label),
        "hurdle": hurdle,
        "hurdle_robust": analyze_hurdle(
            directory,
            label + ":hurdle:robust",
            HURDLE_FEATURE_SETS[hurdle_feature_set],
            ROBUST_TARGET,
        ),
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
    # This is the only geometry/model/feature selection.  Odd-event values do
    # not alter this choice; they motivated only the fixed robust cut target.
    discovery_rows = [
        row for row in rows
        if row["family"] in {"standard_full", "compact", "hurdle"}
    ]
    discovery_selected = max(
        discovery_rows, key=lambda row: row["validation"]["e_keep"]
    )
    if discovery_selected["family"] == "hurdle":
        selected = next(
            row for row in rows
            if row["geometry"] == discovery_selected["geometry"]
            and row["family"] == "hurdle_robust"
            and row["feature_set"] == discovery_selected["feature_set"]
        )
    else:
        selected = discovery_selected
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
        "robust_target": {
            "pim_reject": ROBUST_TARGET[0],
            "pi0_reject": ROBUST_TARGET[1],
            "definition": (
                "fixed after the discovery candidate missed the seed-6302026 "
                "odd background rates and before reading seed 7302026"
            ),
        },
        "split": {
            "fit": "eventID mod 4 == 0",
            "validation": "eventID mod 4 == 2; geometry/method/feature/cut selection",
            "candidate_evaluation": (
                "odd eventID; inspected after discovery selection to set the "
                "robust cut margin, not to revise geometry/model/features"
            ),
        },
        "discovery_selection_rule": (
            "maximum validation electron retention among all declared geometry, "
            "method-family, and feature-set candidates at the exact target"
        ),
        "discovery_selected": discovery_selected,
        "final_operating_rule": (
            "keep the discovery geometry/model/features and retune only its cut "
            "to the fixed robust validation target"
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
