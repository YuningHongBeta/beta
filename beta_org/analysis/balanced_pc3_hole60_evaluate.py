#!/usr/bin/env python3
"""Evaluate T3 Ncluster==1 plus downstream-PC no-hit for the balanced layouts."""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from bgo_features_v2 import load_bgo2
from bgoegg_pc_aperture_evaluate import ensure_features

BASE = Path(__file__).resolve().parents[1]
ROOT_DIR = BASE / "build_bgoegg_frustum" / "output" / "scan"
CACHE_DIR = BASE / "tmp" / "balanced_pc3_hole60_s7232026"
RESULT = BASE / "analysis" / "results" / "bgoegg_v1" / \
    "balanced_pc3_hole60_s7232026.json"
EXTRACTOR = BASE / "analysis" / "bgo_extract_features_v2"

PC_THRESHOLD_MEV = 0.2
SAMPLES = {
    "BGOC": {
        species: ROOT_DIR / f"bgoc_balanced_pc3_h60_s7232026_{species}.root"
        for species in ("e", "pim", "pi0")
    },
    "BGOegg31": {
        species: ROOT_DIR / f"bgoegg31_balanced_pc3_h60_s7232026_{species}.root"
        for species in ("e", "pim", "pi0")
    },
}


def wilson_95(k: int, n: int) -> tuple[float, float]:
    z = 1.959963984540054
    p = k / n
    d = 1.0 + z * z / n
    c = (p + z * z / (2.0 * n)) / d
    h = z * math.sqrt(p * (1.0 - p) / n + z * z / (4.0 * n * n)) / d
    return c - h, c + h


def cluster_hist(values: np.ndarray) -> list[int]:
    return [
        int(np.count_nonzero(values == i)) for i in range(5)
    ] + [int(np.count_nonzero(values >= 5))]


def evaluate(root_file: Path) -> dict:
    cache = CACHE_DIR / f"{root_file.stem}.bgo2"
    ensure_features(root_file.resolve(), cache.resolve(), EXTRACTOR.resolve())
    data, metadata, manifest = load_bgo2(cache, copy=False)
    names = list(manifest["features"])
    col = {name: index for index, name in enumerate(names)}
    ncluster = data[:, col["nCl4"]].astype(np.int64)
    pc_edep = data[:, col["pcUpE_MeV"]]
    bgo_keep = ncluster == 1
    pc_nohit = pc_edep < PC_THRESHOLD_MEV
    selected = bgo_keep & pc_nohit
    n = len(data)
    n_selected = int(np.count_nonzero(selected))
    lo, hi = wilson_95(n_selected, n)
    n_bgo = int(np.count_nonzero(bgo_keep))
    conditional_pc_veto = int(np.count_nonzero(bgo_keep & ~pc_nohit))
    return {
        "root": str(root_file),
        "cache": str(cache),
        "n": n,
        "n_bgo_ncluster1": n_bgo,
        "n_selected": n_selected,
        "selected_fraction": n_selected / n,
        "selected_wilson95": [lo, hi],
        "rejected_fraction": 1.0 - n_selected / n,
        "rejected_wilson95": [1.0 - hi, 1.0 - lo],
        "conditional_pc_veto_count": conditional_pc_veto,
        "conditional_pc_veto_fraction": conditional_pc_veto / n_bgo if n_bgo else None,
        "ncluster_hist_0_1_2_3_4_5plus": cluster_hist(ncluster),
        "pc_edep_quantile_mev": {
            str(q): float(np.quantile(pc_edep, q))
            for q in (0.0, 0.5, 0.9, 0.99, 1.0)
        },
        "threshold_mev_per_bgo_cell": manifest["threshold_MeV"],
        "pc_threshold_mev": PC_THRESHOLD_MEV,
        "geometry_mode": manifest["geometryMode"],
        "n_layer": manifest["nLayer"],
        "n_sector": manifest["nSector"],
    }


def main() -> None:
    missing = [str(path) for samples in SAMPLES.values()
               for path in samples.values() if not path.exists()]
    if missing:
        raise SystemExit("Missing ROOT files:\n" + "\n".join(missing))
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    result = {
        "tag": "balanced_pc3_hole60_s7232026",
        "selection": "T3 nCl4 == 1 and physical-downstream PC Edep < 0.2 MeV",
        "pc": {
            "side": "physical downstream only (legacy beta branch PCUp)",
            "layers": 3,
            "per_layer": "Pb 6 mm + plastic 4 mm",
            "square_hole_mm": 60,
        },
        "samples": {
            geometry: {
                species: evaluate(path)
                for species, path in samples.items()
            }
            for geometry, samples in SAMPLES.items()
        },
    }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(RESULT)
    for geometry, samples in result["samples"].items():
        e = samples["e"]
        pim = samples["pim"]
        pi0 = samples["pi0"]
        print(
            f"{geometry}: e keep={100*e['selected_fraction']:.3f}% "
            f"pim reject={100*pim['rejected_fraction']:.3f}% "
            f"pi0 reject={100*pi0['rejected_fraction']:.3f}%"
        )


if __name__ == "__main__":
    main()
