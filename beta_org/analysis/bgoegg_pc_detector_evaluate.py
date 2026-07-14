#!/usr/bin/env python3
"""Evaluate each beta-decay detector and a stand-alone BGOegg photon veto.

The Photon Counter (PC) is never used as a classifier feature here.  A PC veto
is a fixed cut on the summed plastic energy after fixed thesis-like detector
selections have been evaluated independently.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import yaml

from bgo_features_v2 import load_bgo2


PROJECT = Path(__file__).resolve().parent.parent
DEFAULT_INPUT = PROJECT / "tmp" / "bgoegg_pc_design_s6302026"
DEFAULT_MANIFEST = (
    PROJECT / "runmanager" / "examples" / "bgoegg_pc_design_s6302026.yml"
)
DEFAULT_OUTPUT = (
    PROJECT / "analysis" / "results" / "bgoegg_v1" /
    "bgoegg_pc_detector_evaluation_s6302026.json"
)
GEOMETRIES = (
    "up4_o12_pb1", "up8_o11p5_pb1", "up8_o12_pb1",
    "up8_o13_pb1", "up12_o12_pb1", "up8_o12_pb2",
    "two8_o12_pb1", "down8_leak_pb1",
)
PC_THRESHOLDS_MEV = (0.0, 0.1, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0)
SPECIES = ("electron", "pim", "pi0")
FILES = {"electron": "e.bgo2", "pim": "pim.bgo2", "pi0": "pi0.bgo2"}


def wilson(p: float, n: int) -> list[float] | None:
    if n <= 0:
        return None
    z = 1.959963984540054
    denominator = 1.0 + z * z / n
    center = (p + z * z / (2.0 * n)) / denominator
    half = z * math.sqrt(
        p * (1.0 - p) / n + z * z / (4.0 * n * n)
    ) / denominator
    return [float(center - half), float(center + half)]


def fraction(mask: np.ndarray) -> dict:
    value = float(np.mean(mask))
    return {
        "value": value,
        "numerator": int(np.count_nonzero(mask)),
        "denominator": int(len(mask)),
        "wilson_95ci": wilson(value, len(mask)),
    }


def conditional_fraction(mask: np.ndarray, condition: np.ndarray) -> dict:
    denominator = int(np.count_nonzero(condition))
    numerator = int(np.count_nonzero(mask & condition))
    if denominator == 0:
        return {
            "value": None, "numerator": numerator, "denominator": 0,
            "wilson_95ci": None,
        }
    value = numerator / denominator
    return {
        "value": float(value), "numerator": numerator,
        "denominator": denominator, "wilson_95ci": wilson(value, denominator),
    }


def selection_metrics(masks: dict[str, np.ndarray]) -> dict:
    return {
        "electron_keep": fraction(masks["electron"]),
        "pim_reject": fraction(~masks["pim"]),
        "pi0_reject": fraction(~masks["pi0"]),
    }


def conditional_selection_metrics(
    masks: dict[str, np.ndarray], conditions: dict[str, np.ndarray]
) -> dict:
    return {
        "electron_keep": conditional_fraction(
            masks["electron"], conditions["electron"]),
        "pim_reject": conditional_fraction(
            ~masks["pim"], conditions["pim"]),
        "pi0_reject": conditional_fraction(
            ~masks["pi0"], conditions["pi0"]),
    }


def veto_metrics(hit: dict[str, np.ndarray], condition=None) -> dict:
    if condition is None:
        return {
            "electron_false_veto": fraction(hit["electron"]),
            "pim_veto": fraction(hit["pim"]),
            "pi0_veto": fraction(hit["pi0"]),
        }
    return {
        "electron_false_veto": conditional_fraction(
            hit["electron"], condition["electron"]),
        "pim_veto": conditional_fraction(hit["pim"], condition["pim"]),
        "pi0_veto": conditional_fraction(hit["pi0"], condition["pi0"]),
    }


def load_geometry(directory: Path) -> tuple[dict, list[str], dict]:
    samples, names, manifests = {}, None, {}
    metadata = None
    for species in SPECIES:
        data, meta, manifest = load_bgo2(directory / FILES[species], copy=False)
        these_names = list(manifest["features"])
        if names is None:
            names = these_names
            metadata = meta
        elif names != these_names:
            raise ValueError(f"feature mismatch in {directory}")
        if not math.isclose(float(meta["threshold_MeV"]), 3.0):
            raise ValueError(f"BGO threshold is not 3 MeV in {directory}")
        samples[species] = data
        manifests[species] = manifest
    assert names is not None and metadata is not None
    return samples, names, {"metadata": metadata, "manifests": manifests}


def detector_masks(samples: dict, names: list[str]) -> dict[str, dict]:
    col = {name: index for index, name in enumerate(names)}
    required = {
        "nCl4", "thTotalE_MeV", "thGeomPath_mm",
        "thMatchedDEdx_MeV_per_mm", "tlcAnySegHitGe1_eff05",
        "absDeltaZValid", "absDeltaZLt90",
    }
    missing = sorted(required - set(col))
    if missing:
        raise ValueError(f"required detector features missing: {missing}")
    output = {name: {} for name in (
        "bgo_ncluster1", "th_total_dedx", "tlc_hit",
        "delta_z", "pre_pc_sequence",
    )}
    for species, data in samples.items():
        bgo = data[:, col["nCl4"]] == 1
        th_valid = data[:, col["thGeomPath_mm"]] > 0.0
        th = (
            th_valid
            & (data[:, col["thTotalE_MeV"]] >= 0.7)
            & (data[:, col["thTotalE_MeV"]] <= 3.0)
            & (data[:, col["thMatchedDEdx_MeV_per_mm"]] >= 0.12)
            & (data[:, col["thMatchedDEdx_MeV_per_mm"]] <= 0.25)
        )
        tlc = data[:, col["tlcAnySegHitGe1_eff05"]] > 0.5
        delta_z = (
            (data[:, col["absDeltaZValid"]] > 0.5)
            & (data[:, col["absDeltaZLt90"]] > 0.5)
        )
        output["bgo_ncluster1"][species] = bgo
        output["th_total_dedx"][species] = th
        output["tlc_hit"][species] = tlc
        output["delta_z"][species] = delta_z
        output["pre_pc_sequence"][species] = bgo & th & tlc & delta_z
    return output


def tlc_scan(samples: dict, names: list[str]) -> list[dict]:
    col = {name: index for index, name in enumerate(names)}
    rows = []
    for efficiency in (1, 2, 5, 10):
        for npe in (1, 2, 3):
            feature = f"tlcAnySegHitGe{npe}_eff{efficiency:02d}"
            if feature not in col:
                continue
            masks = {
                species: data[:, col[feature]] > 0.5
                for species, data in samples.items()
            }
            rows.append({
                "collection_times_pde_percent": efficiency,
                "segment_npe_threshold": npe,
                "selection": "electron candidate requires any TLC segment hit",
                "metrics_all_generated": selection_metrics(masks),
            })
    return rows


def cutflow(detectors: dict[str, dict[str, np.ndarray]]) -> list[dict]:
    cumulative = {species: np.ones_like(
        detectors["bgo_ncluster1"][species], dtype=bool) for species in SPECIES}
    rows = []
    for stage in ("bgo_ncluster1", "th_total_dedx", "tlc_hit", "delta_z"):
        before = {species: cumulative[species].copy() for species in SPECIES}
        for species in SPECIES:
            cumulative[species] &= detectors[stage][species]
        rows.append({
            "stage": stage,
            "selection_metrics_all_generated": selection_metrics(cumulative),
            "conditional_pass_given_previous": {
                species: conditional_fraction(cumulative[species], before[species])
                for species in SPECIES
            },
        })
    return rows


def pc_hit(data: np.ndarray, index: int, threshold: float) -> np.ndarray:
    if threshold == 0.0:
        return data[:, index] > 0.0
    return data[:, index] >= threshold


def evaluate_geometry(directory: Path, name: str, design: dict) -> dict:
    samples, names, provenance = load_geometry(directory)
    col = {feature: index for index, feature in enumerate(names)}
    detectors = detector_masks(samples, names)
    delta_z_valid = {
        species: data[:, col["absDeltaZValid"]] > 0.5
        for species, data in samples.items()
    }
    response = []
    for threshold in PC_THRESHOLDS_MEV:
        total_hit = {
            species: pc_hit(data, col["pcSumE_MeV"], threshold)
            for species, data in samples.items()
        }
        up_hit = {
            species: pc_hit(data, col["pcUpE_MeV"], threshold)
            for species, data in samples.items()
        }
        down_hit = {
            species: pc_hit(data, col["pcDownE_MeV"], threshold)
            for species, data in samples.items()
        }
        response.append({
            "pc_plastic_threshold_MeV": threshold,
            "hit_rule": "EdepPC>0" if threshold == 0.0 else (
                f"EdepPC >= {threshold:g} MeV"),
            "all_generated_pc_only": veto_metrics(total_hit),
            "given_bgo_ncluster1": veto_metrics(
                total_hit, detectors["bgo_ncluster1"]),
            "given_full_pre_pc_sequence": veto_metrics(
                total_hit, detectors["pre_pc_sequence"]),
            "side_response_all_generated": {
                "upstream": veto_metrics(up_hit),
                "downstream": veto_metrics(down_hit),
            },
        })
    return {
        "design": design,
        "source_root_files": {
            species: provenance["manifests"][species]["input"]
            for species in SPECIES
        },
        "bgo2_metadata": provenance["metadata"],
        "fixed_detector_definitions": {
            "bgo": "Nclus4 == 1; per-segment BGO threshold fixed at 3 MeV",
            "th": (
                "require valid target-to-leading-BGO geometric path, "
                "0.7<=total dE<=3.0 MeV and 0.12<=dE/dx<=0.25 MeV/mm"
            ),
            "tlc": (
                "any of 30 segments Npe>=1 at collection*PDE=5%; "
                "provisional optical-response point"
            ),
            "delta_z": "require valid TH z readout and abs(zReco-zPred)<90 mm",
        },
        "independent_detector_performance_all_generated": {
            stage: selection_metrics(detectors[stage])
            for stage in ("bgo_ncluster1", "th_total_dedx", "tlc_hit", "delta_z")
        },
        "delta_z_performance_given_valid_readout": (
            conditional_selection_metrics(detectors["delta_z"], delta_z_valid)
        ),
        "tlc_response_scan_all_generated": tlc_scan(samples, names),
        "fixed_cutflow": cutflow(detectors),
        "pre_pc_selection_all_generated": selection_metrics(
            detectors["pre_pc_sequence"]),
        "pc_standalone_veto": response,
    }


def point_at(result: dict, threshold: float) -> dict:
    return next(row for row in result["pc_standalone_veto"]
                if row["pc_plastic_threshold_MeV"] == threshold)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    manifest = yaml.safe_load(args.manifest.read_text(encoding="utf-8"))
    designs = {item["name"]: item for item in manifest["geometries"]}
    geometries = {
        name: evaluate_geometry(args.input / name, name, designs[name])
        for name in GEOMETRIES
    }

    # A stable screen with useful statistics: hard PC-only veto after the fixed
    # BGO Nclus=1 selection.  The much smaller full pre-PC denominator is still
    # reported and remains the thesis-comparable primary metric.
    screen_threshold = 0.5
    screen_rows = []
    for name, result in geometries.items():
        point = point_at(result, screen_threshold)
        screen_rows.append({
            "geometry": name,
            "pc_plastic_threshold_MeV": screen_threshold,
            "given_bgo_ncluster1": point["given_bgo_ncluster1"],
            "given_full_pre_pc_sequence": point["given_full_pre_pc_sequence"],
            "all_generated_pc_only": point["all_generated_pc_only"],
        })
    selected = max(
        screen_rows,
        key=lambda row: (
            row["given_bgo_ncluster1"]["pi0_veto"]["value"] or -1.0,
            -(row["given_bgo_ncluster1"]["electron_false_veto"]["value"] or 0.0),
        ),
    )
    result = {
        "status": "complete",
        "analysis": "bgoegg_pc_detector_evaluation",
        "selection_policy": (
            "PC is a hard stand-alone plastic-energy veto and is never a BGO "
            "classifier input. The 0.5 MeV screening point ranks designs by pi0 "
            "veto conditional on the fixed BGO Nclus=1 selection; the full "
            "thesis-like pre-PC conditional value is reported separately and "
            "is not optimized because its 100k denominator can be small."
        ),
        "bgo_threshold_MeV": 3.0,
        "thesis_benchmarks": {
            "hong_table_2_2_remaining_fraction": {
                "bgo_pi0": 0.03, "bgo_pim": 0.07,
                "th_pi0": 0.10, "th_pim": 0.001,
                "tlc_pi0": 0.10, "tlc_pim": 0.045,
                "delta_z_pi0": 0.80, "pc_pi0_given_pre_pc": 0.30,
            },
            "kamada_pc": {
                "ideal_geometry_veto_given_remaining_pi0": 0.70,
                "realistic_100MeV_gamma_detection_max": 0.64,
                "caveat": (
                    "Kamada section 3.4 and summary disagree on the Pb thickness "
                    "assigned to the two-layer 64% point"
                ),
            },
        },
        "pc_response_scope": {
            "measured_here": (
                "effective PC-only veto including geometric acceptance and plastic "
                "energy deposit"
            ),
            "not_separated_in_current_root": (
                "gamma intersection with the first PC layer versus intrinsic "
                "detection conditional on intersection"
            ),
            "real_background_occupancy": "not simulated",
        },
        "screen_threshold_MeV": screen_threshold,
        "selected_screening_design": selected,
        "screening_table": screen_rows,
        "geometries": geometries,
        "design_manifest": str(args.manifest.resolve()),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "selected_screening_design": selected,
        "full_pre_pc_pi0_denominators": {
            name: point_at(value, screen_threshold)[
                "given_full_pre_pc_sequence"]["pi0_veto"]["denominator"]
            for name, value in geometries.items()
        },
    }, indent=2))
    print(args.output)


if __name__ == "__main__":
    main()
