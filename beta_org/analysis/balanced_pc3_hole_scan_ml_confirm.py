#!/usr/bin/env python3
"""Apply frozen hole-specific BR-weighted BDT and fixed cut to confirmation BGO2."""

from __future__ import annotations

import argparse
from array import array
import json
import math
from pathlib import Path

import numpy as np
import ROOT

from bgo_features_v2 import load_bgo2


PI0_TO_BETA = 3.33333333333333
PIM_TO_BETA = 0.0375
SPECIES = ("e", "pim", "pi0")
PREFIX = {"BGOC": "bgoc", "BGOegg31": "bgoegg31"}


def wilson_upper(selected: int, denominator: int) -> float:
    z = 1.959963984540054
    p = selected / denominator
    scale = 1.0 + z * z / denominator
    centre = (p + z * z / (2.0 * denominator)) / scale
    half = z * math.sqrt(
        p * (1.0 - p) / denominator + z * z / (4.0 * denominator**2)
    ) / scale
    return centre + half


def summarize(selected: dict[str, int], denominator: int) -> dict:
    pi0_survive = selected["pi0"] / denominator
    pim_survive = selected["pim"] / denominator
    return {
        "denominator_per_species": denominator,
        "electron_keep_count": selected["e"],
        "electron_keep": selected["e"] / denominator,
        "pim_reject_count": denominator - selected["pim"],
        "pim_reject": 1.0 - pim_survive,
        "pi0_reject_count": denominator - selected["pi0"],
        "pi0_reject": 1.0 - pi0_survive,
        "pi_total_to_beta": PI0_TO_BETA * pi0_survive + PIM_TO_BETA * pim_survive,
        "pi_total_to_beta_mc_stat_upper95": (
            PI0_TO_BETA * wilson_upper(selected["pi0"], denominator)
            + PIM_TO_BETA * wilson_upper(selected["pim"], denominator)
        ),
    }


def load(path: Path) -> tuple[np.ndarray, dict[str, int]]:
    values, _, manifest = load_bgo2(path, copy=False)
    if len(values) != 100_000:
        raise ValueError(f"{path}: expected 100000 rows, found {len(values)}")
    return values, {name: index for index, name in enumerate(manifest["features"])}


def sample_path(input_dir: Path, geometry: str, sample_tag: str,
                species: str) -> Path:
    return input_dir / f"{PREFIX[geometry]}_{sample_tag}_{species}.bgo2"


def fixed_cut(input_dir: Path, geometry: str, sample_tag: str) -> dict:
    selected = {}
    for species in SPECIES:
        data, column = load(sample_path(input_dir, geometry, sample_tag, species))
        mask = (data[:, column["nCl4"]] == 1) & (data[:, column["pcUpE_MeV"]] < 0.2)
        selected[species] = int(np.count_nonzero(mask))
    return summarize(selected, 100_000)


def weighted_bdt(input_dir: Path, geometry: str, sample_tag: str,
                 development: dict) -> tuple[dict, dict, str, dict]:
    model = development["geometries"][geometry]["models"]["bgo_pc_brweighted"]
    thresholds = {
        "central": float(model["br4_central"]["validation"]["threshold"]),
        "guarded": float(model["br4_guarded"]["validation"]["threshold"]),
    }
    weight_path = Path(model["weight_file"])
    reader = ROOT.TMVA.Reader("!Color:!Silent")
    buffers = {name: array("f", [0.0]) for name in model["features"]}
    for name, value in buffers.items():
        reader.AddVariable(name, value)
    reader.BookMVA(development["model"]["method"], str(weight_path))
    selected = {name: {} for name in thresholds}
    score_edges = np.linspace(-1.0, 1.0, 101)
    score_histograms = {}
    for species in SPECIES:
        data, column = load(sample_path(input_dir, geometry, sample_tag, species))
        scores = np.empty(len(data), dtype=np.float64)
        for index, row in enumerate(data):
            for name, value in buffers.items():
                value[0] = float(row[column[name]])
            score = reader.EvaluateMVA(development["model"]["method"])
            scores[index] = score
        for name, threshold in thresholds.items():
            selected[name][species] = int(np.count_nonzero(scores > threshold))
        histogram, _ = np.histogram(scores, bins=score_edges)
        score_histograms[species] = {
            "counts": histogram.tolist(),
            "underflow": int(np.count_nonzero(scores < score_edges[0])),
            "overflow": int(np.count_nonzero(scores > score_edges[-1])),
        }
    return ({name: summarize(counts, 100_000)
             for name, counts in selected.items()},
            thresholds, str(weight_path), {
        "edges": score_edges.tolist(),
        "species": score_histograms,
    })


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir", type=Path)
    parser.add_argument("development_json", type=Path)
    parser.add_argument("confirmation_sample_tag")
    parser.add_argument("output_json", type=Path)
    args = parser.parse_args()
    development = json.loads(args.development_json.read_text())
    results = {}
    for geometry in PREFIX:
        fixed = fixed_cut(args.input_dir, geometry, args.confirmation_sample_tag)
        bdt, thresholds, weight, score_histogram = weighted_bdt(
            args.input_dir, geometry, args.confirmation_sample_tag, development
        )
        results[geometry] = {
            "fixed_t3_pc": fixed,
            "brweighted_bdt_central": bdt["central"],
            "brweighted_bdt_guarded": bdt["guarded"],
            "frozen_thresholds": thresholds,
            "weight_file": weight,
            "confirmation_score_histogram": score_histogram,
        }
        for name, metrics in (("fixed", fixed),
                              ("bdt-central", bdt["central"]),
                              ("bdt-guarded", bdt["guarded"])):
            print(
                f"{geometry} {name}: e={100*metrics['electron_keep']:.3f}% "
                f"pimR={100*metrics['pim_reject']:.3f}% "
                f"pi0R={100*metrics['pi0_reject']:.3f}% "
                f"Rbg={100*metrics['pi_total_to_beta']:.3f}% "
                f"upper95={100*metrics['pi_total_to_beta_mc_stat_upper95']:.3f}%"
            )
    output = {
        "status": "independent-seed frozen-cut confirmation",
        "development_json": str(args.development_json),
        "development_seed": 7232026,
        "confirmation_seed": 7242026,
        "confirmation_sample_tag": args.confirmation_sample_tag,
        "events_per_species": 100_000,
        "results": results,
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(output, indent=2) + "\n")


if __name__ == "__main__":
    main()
