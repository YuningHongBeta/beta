#!/usr/bin/env python3
"""Plot thesis-aware detector and stand-alone PC performance."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages


def read(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def point(geometry: dict, threshold: float = 0.5) -> dict:
    return next(row for row in geometry["operating_points"]
                if row["pc_plastic_threshold_MeV"] == threshold)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--detector", type=Path, required=True)
    parser.add_argument("--aperture", type=Path, required=True)
    parser.add_argument("--wide", type=Path, required=True)
    parser.add_argument("--annulus", type=Path, required=True)
    parser.add_argument("--confirm8", type=Path, required=True)
    parser.add_argument("--confirm12", type=Path, required=True)
    parser.add_argument("--pdf", type=Path, required=True)
    parser.add_argument("--png", type=Path, required=True)
    args = parser.parse_args()

    detector = read(args.detector)
    aperture = read(args.aperture)
    wide = read(args.wide)
    annulus = read(args.annulus)
    confirms = [read(args.confirm8), read(args.confirm12)]

    baseline = detector["geometries"]["up8_o12_pb2"]
    independent = baseline["independent_detector_performance_all_generated"]
    dz_valid = baseline["delta_z_performance_given_valid_readout"]
    current_pi0 = [
        independent["bgo_ncluster1"]["pi0_reject"]["value"],
        independent["th_total_dedx"]["pi0_reject"]["value"],
        independent["tlc_hit"]["pi0_reject"]["value"],
        dz_valid["pi0_reject"]["value"],
    ]
    current_pim = [
        independent["bgo_ncluster1"]["pim_reject"]["value"],
        independent["th_total_dedx"]["pim_reject"]["value"],
        independent["tlc_hit"]["pim_reject"]["value"],
        dz_valid["pim_reject"]["value"],
    ]
    thesis_pi0 = [0.97, 0.90, 0.90, 0.20]
    thesis_pim = [0.93, 0.999, 0.955, np.nan]

    fig, axes = plt.subplots(2, 2, figsize=(14.0, 9.3), constrained_layout=True)
    fig.patch.set_facecolor("white")

    ax = axes[0, 0]
    x = np.arange(4)
    width = 0.37
    ax.bar(x - width / 2, np.array(thesis_pi0) * 100, width,
           label="Hong thesis reference", color="#BAB0AC")
    ax.bar(x + width / 2, np.array(current_pi0) * 100, width,
           label="exact 31-ring BGOegg MC", color="#4E79A7")
    for index, value in enumerate(current_pi0):
        ax.text(index + width / 2, value * 100 - 4.5, f"{value*100:.1f}",
                ha="center", va="top", fontsize=8, color="white")
    ax.set_xticks(x, ["BGO", "TH", "TLC", r"$\Delta z$"])
    ax.set_ylim(0, 108)
    ax.set_ylabel(r"$\pi^0$ rejection (%)")
    ax.set_title("A. Individual detector response (not a joint classifier)")
    ax.grid(axis="y", alpha=0.2)
    ax.legend(fontsize=8, loc="lower left")
    ax.text(0.01, 0.98,
            "BGO/TH/TLC: all generated; Δz: valid-readout denominator\n"
            "BGO per-segment threshold = 3 MeV",
            transform=ax.transAxes, va="top", fontsize=8,
            bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "none"})

    ax = axes[0, 1]
    labels = ["annulus\n12×Pb3", "holes\n4×Pb3",
              "wide\n4×Pb3", "wide\n8×Pb2", "wide\n12×Pb2"]
    sources = [
        annulus["geometries"]["two_annulus_12_pb3"],
        aperture["geometries"]["two_full_4_pb3"],
        wide["geometries"]["two_wide_4_pb3"],
        wide["geometries"]["two_wide_8_pb2"],
        wide["geometries"]["two_wide_12_pb2"],
    ]
    intercept, response, effective = [], [], []
    for source in sources:
        q = point(source)["given_bgo_ncluster1"]["pi0"]
        intercept.append(q["gamma_intercept"]["value"] * 100)
        response.append(q["intrinsic_hit_given_gamma_intercept"]["value"] * 100)
        effective.append(q["effective_pc_hit"]["value"] * 100)
    x = np.arange(len(labels)); width = 0.25
    ax.bar(x - width, intercept, width, label="gamma intercept", color="#59A14F")
    ax.bar(x, response, width, label="hit | intercept", color="#F28E2B")
    ax.bar(x + width, effective, width, label="effective hard veto", color="#E15759")
    ax.axhline(70, color="black", linestyle="--", linewidth=1.1,
               label="thesis PC target 70%")
    ax.set_xticks(x, labels)
    ax.set_ylim(0, 108)
    ax.set_ylabel(r"rate given BGO $N_{cluster}=1$ (%)")
    ax.set_title("B. PC acceptance and response decomposition (0.5 MeV)")
    ax.grid(axis="y", alpha=0.2)
    ax.legend(fontsize=7.5, ncol=2, loc="lower right")

    ax = axes[1, 0]
    names = ["wide 8×(Pb2+plastic5)", "wide 12×(Pb2+plastic5)"]
    values, lows, highs = [], [], []
    for result in confirms:
        metric = result["fixed_point"]["conditions"][
            "given_full_pre_pc"]["effective_pc_hit"]
        values.append(metric["value"] * 100)
        lows.append(metric["wilson_95ci"][0] * 100)
        highs.append(metric["wilson_95ci"][1] * 100)
    yerr = np.array([np.array(values) - np.array(lows),
                     np.array(highs) - np.array(values)])
    bars = ax.bar(np.arange(2), values, color=["#4E79A7", "#E15759"],
                  yerr=yerr, capsize=5)
    ax.axhline(70, color="black", linestyle="--", linewidth=1.1)
    ax.set_xticks(np.arange(2), names)
    ax.set_ylim(60, 95)
    ax.set_ylabel(r"PC-only $\pi^0$ veto after all pre-PC cuts (%)")
    ax.set_title("C. Independent-seed 2M-event confirmation")
    ax.grid(axis="y", alpha=0.2)
    for bar, value, result in zip(bars, values, confirms):
        n = result["fixed_point"]["conditions"][
            "given_full_pre_pc"]["effective_pc_hit"]
        ax.text(bar.get_x() + bar.get_width()/2, value + 1.0,
                f"{value:.1f}% ({n['numerator']}/{n['denominator']})",
                ha="center", fontsize=8.5)
    ax.text(0.01, 0.03,
            "95% Wilson intervals; wide disks block charged-particle corridors\n"
            "Hong PC 70% is conditional on the remaining pre-PC π0 sample",
            transform=ax.transAxes, fontsize=8)

    ax = axes[1, 1]
    for result, label, color in zip(
            confirms, ("wide 8 layers", "wide 12 layers"),
            ("#4E79A7", "#E15759")):
        thresholds, rates = [], []
        for row in result["operating_points"]:
            thresholds.append(row["pc_plastic_threshold_MeV"])
            rates.append(row["conditions"]["given_full_pre_pc"]
                         ["effective_pc_hit"]["value"] * 100)
        ax.plot(thresholds, rates, marker="o", label=label, color=color)
    ax.axhline(70, color="black", linestyle="--", linewidth=1.1,
               label="thesis target")
    ax.axvline(0.5, color="#777777", linestyle=":", linewidth=1.0,
               label="fixed PC point")
    ax.set_xlim(-0.2, 10.2); ax.set_ylim(0, 95)
    ax.set_xlabel("PC plastic energy threshold (MeV)")
    ax.set_ylabel(r"PC-only $\pi^0$ veto after all pre-PC cuts (%)")
    ax.set_title("D. PC threshold dependence; BGO remains fixed at 3 MeV")
    ax.grid(alpha=0.2)
    ax.legend(fontsize=8, loc="lower left")

    wide8_e = point(wide["geometries"]["two_wide_8_pb2"])[
        "given_bgo_ncluster1"]["e"]["effective_pc_hit"]["value"] * 100
    wide12_e = point(wide["geometries"]["two_wide_12_pb2"])[
        "given_bgo_ncluster1"]["e"]["effective_pc_hit"]["value"] * 100
    fig.suptitle(
        "Exact 31-ring BGOegg detector audit and stand-alone Photon Counter optimization\n"
        f"BGO threshold 3 MeV; PC hard-hit veto; electron false veto after BGO: "
        f"8 layers {wide8_e:.3f}%, 12 layers {wide12_e:.3f}%",
        fontsize=14,
    )
    args.pdf.parent.mkdir(parents=True, exist_ok=True)
    args.png.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(args.pdf) as pdf:
        pdf.savefig(fig, bbox_inches="tight", facecolor="white")
    fig.savefig(args.png, dpi=180, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(args.pdf)
    print(args.png)


if __name__ == "__main__":
    main()
