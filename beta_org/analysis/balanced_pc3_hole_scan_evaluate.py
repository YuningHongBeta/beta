#!/usr/bin/env python3
"""Evaluate square/circular PC aperture scans with fixed T3+PC selection."""

from __future__ import annotations

import json
import math
import os
from pathlib import Path

import numpy as np

from bgo_features_v2 import load_bgo2
from bgoegg_pc_aperture_evaluate import ensure_features


BASE = Path(__file__).resolve().parents[1]
ROOT_DIR = BASE / "build_bgoegg_frustum/output/scan"
CACHE_DIR = Path(os.environ.get(
    "BETA_HOLE_SCAN_CACHE",
    "/gpfs/group/had/sks/Users/yhong/rootfile/beta/"
    "balanced_pc3_hole_scan_s7232026",
))
RESULT = BASE / "analysis/results/bgoegg_v1/balanced_pc3_hole_scan_20260722.json"
EXTRACTOR = BASE / "analysis/bgo_extract_features_v2"

SPECIES = ("e", "pim", "pi0")
GEOMETRIES = {"BGOC": "bgoc", "BGOegg31": "bgoegg31"}
PI0_TO_BETA = 3.33333333333333
PIM_TO_BETA = 0.0375
PC_THRESHOLD_MEV = 0.2
POSITIVE_PC_EDGES_MEV = np.geomspace(1.0e-4, 50.0, 81)


def wilson_95(k: int, n: int) -> tuple[float, float]:
    z = 1.959963984540054
    p = k / n
    d = 1.0 + z * z / n
    c = (p + z * z / (2.0 * n)) / d
    h = z * math.sqrt(p * (1.0 - p) / n + z * z / (4.0 * n * n)) / d
    return c - h, c + h


def root_path(shape: str, diameter_or_side_mm: int, seed: int,
              geometry_prefix: str, species: str) -> Path:
    hole_tag = f"hc{diameter_or_side_mm}" if shape == "circle" else f"h{diameter_or_side_mm}"
    if shape == "square" and diameter_or_side_mm == 60:
        confirmation = "_confirm" if seed == 7242026 else ""
        name = f"{geometry_prefix}_balanced_pc3_h60{confirmation}_s{seed}_{species}.root"
    else:
        name = f"{geometry_prefix}_balanced_pc3_{hole_tag}_s{seed}_{species}.root"
    return ROOT_DIR / name


def evaluate(path: Path) -> dict:
    cache = CACHE_DIR / f"{path.stem}.bgo2"
    ensure_features(path.resolve(), cache.resolve(), EXTRACTOR.resolve())
    data, _, manifest = load_bgo2(cache, copy=False)
    if len(data) != 100_000:
        raise ValueError(f"{path}: expected 100000 rows, found {len(data)}")
    column = {name: index for index, name in enumerate(manifest["features"])}
    ncluster = data[:, column["nCl4"]].astype(np.int64)
    pc_edep = data[:, column["pcUpE_MeV"]]
    bgo_keep = ncluster == 1
    selected = bgo_keep & (pc_edep < PC_THRESHOLD_MEV)
    positive = pc_edep[bgo_keep & (pc_edep > 0)]
    positive_hist, _ = np.histogram(positive, bins=POSITIVE_PC_EDGES_MEV)
    n_selected = int(np.count_nonzero(selected))
    lo, hi = wilson_95(n_selected, len(data))
    return {
        "root": str(path),
        "cache": str(cache),
        "n": len(data),
        "n_bgo_ncluster1": int(np.count_nonzero(bgo_keep)),
        "n_selected": n_selected,
        "selected_fraction": n_selected / len(data),
        "selected_wilson95": [lo, hi],
        "rejected_fraction": 1.0 - n_selected / len(data),
        "ncluster_hist_0_1_2_3_4_5plus": [
            *[int(np.count_nonzero(ncluster == value)) for value in range(5)],
            int(np.count_nonzero(ncluster >= 5)),
        ],
        "pc_after_ncluster1": {
            "zero_count": int(np.count_nonzero(bgo_keep & (pc_edep == 0))),
            "positive_hist_edges_mev": POSITIVE_PC_EDGES_MEV.tolist(),
            "positive_hist_counts": positive_hist.tolist(),
            "positive_overflow_count": int(np.count_nonzero(positive > POSITIVE_PC_EDGES_MEV[-1])),
        },
    }


def summarize_species(samples: dict) -> dict:
    e = samples["e"]
    pim = samples["pim"]
    pi0 = samples["pi0"]
    total = PI0_TO_BETA * pi0["selected_fraction"] + PIM_TO_BETA * pim["selected_fraction"]
    return {
        "electron_keep": e["selected_fraction"],
        "pim_reject": pim["rejected_fraction"],
        "pi0_reject": pi0["rejected_fraction"],
        "pi_total_to_beta": total,
    }


def main() -> None:
    configurations = [
        ("square", 60, seed) for seed in (7232026, 7242026)
    ] + [
        ("square", size, 7232026) for size in (100, 150)
    ] + [
        ("circle", size, seed)
        for seed in (7232026, 7242026) for size in (60, 100, 150)
    ]
    missing = [
        str(root_path(shape, size, seed, prefix, species))
        for shape, size, seed in configurations
        for prefix in GEOMETRIES.values() for species in SPECIES
        if not root_path(shape, size, seed, prefix, species).exists()
    ]
    if missing:
        raise SystemExit("Missing ROOT files:\n" + "\n".join(missing))
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    scans = {}
    for shape, size, seed in configurations:
        key = f"{shape}_{size}mm_s{seed}"
        scans[key] = {
            "shape": shape,
            "diameter_or_side_mm": size,
            "seed": seed,
            "geometries": {},
        }
        for geometry, prefix in GEOMETRIES.items():
            samples = {
                species: evaluate(root_path(shape, size, seed, prefix, species))
                for species in SPECIES
            }
            scans[key]["geometries"][geometry] = {
                "samples": samples,
                "summary": summarize_species(samples),
            }
            summary = scans[key]["geometries"][geometry]["summary"]
            print(
                f"{key} {geometry}: e={100*summary['electron_keep']:.3f}% "
                f"pimR={100*summary['pim_reject']:.3f}% "
                f"pi0R={100*summary['pi0_reject']:.3f}% "
                f"Rbg={100*summary['pi_total_to_beta']:.3f}%"
            )
    result = {
        "tag": "balanced_pc3_hole_scan_20260722",
        "selection": "T3 nCl4 == 1 and physical-downstream PC Edep < 0.2 MeV",
        "branching_requirement": {
            "formula": "3.33333333333333*pi0_survival + 0.0375*pim_survival",
            "maximum": 0.04,
        },
        "pc": {"side": "physical downstream only", "layers": 3,
               "per_layer": "Pb 6 mm + plastic 4 mm"},
        "scans": scans,
    }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(RESULT)


if __name__ == "__main__":
    main()
