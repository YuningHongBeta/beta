#!/usr/bin/env python3
"""Draw the exact 31-layer BGOegg and the selected photon-counter endcap."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Polygon

from plot_bgo_geometry_comparison import CRYSTAL_LENGTH_CM, egg31_rings


Z_OFFSET_CM = -10.0


def egg_polygon(
    theta_low: float, theta_high: float, front: float, x_sign: float,
) -> list[tuple[float, float]]:
    points = []
    for radius, theta in (
        (front, theta_low),
        (front + CRYSTAL_LENGTH_CM, theta_low),
        (front + CRYSTAL_LENGTH_CM, theta_high),
        (front, theta_high),
    ):
        angle = math.radians(theta)
        points.append((
            x_sign * radius * math.sin(angle),
            Z_OFFSET_CM + radius * math.cos(angle),
        ))
    return points


def endcap_polygon(
    z_near: float, z_far: float, theta_inner: float, theta_outer: float,
    z_sign: float, x_sign: float,
) -> list[tuple[float, float]]:
    inner = math.tan(math.radians(theta_inner))
    outer = math.tan(math.radians(theta_outer))
    return [
        (x_sign * z_near * inner, z_sign * z_near),
        (x_sign * z_far * inner, z_sign * z_far),
        (x_sign * z_far * outer, z_sign * z_far),
        (x_sign * z_near * outer, z_sign * z_near),
    ]


def draw_endcap(axis: plt.Axes, design: dict, direction: str) -> None:
    if design["photon_counter"] not in {direction, "two_sided"}:
        return
    if direction == "upstream":
        inner = design["pc_up_theta_inner_deg"]
        outer = design["pc_up_theta_outer_deg"]
        z_sign = -1.0
    else:
        inner = design["pc_down_theta_inner_deg"]
        outer = design["pc_down_theta_outer_deg"]
        z_sign = 1.0
    z_front = design["pc_z_front_cm"]
    pb_cm = 0.1 * design["pc_pb_thickness_mm"]
    scinti_cm = 0.1 * design["pc_scinti_thickness_mm"]
    layer_cm = pb_cm + scinti_cm
    for layer in range(design["pc_n_layers"]):
        pb_near = z_front + layer * layer_cm
        pb_far = pb_near + pb_cm
        scinti_far = pb_far + scinti_cm
        for x_sign in (-1.0, 1.0):
            axis.add_patch(Polygon(
                endcap_polygon(
                    pb_near, pb_far, inner, outer, z_sign, x_sign
                ),
                closed=True, facecolor="#7F7F7F", edgecolor="none",
            ))
            axis.add_patch(Polygon(
                endcap_polygon(
                    pb_far, scinti_far, inner, outer, z_sign, x_sign
                ),
                closed=True, facecolor="#F2CF5B", edgecolor="none",
            ))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--result", type=Path, required=True)
    parser.add_argument("--geometry")
    parser.add_argument("--pdf", type=Path, required=True)
    parser.add_argument("--png", type=Path, required=True)
    args = parser.parse_args()
    result = json.loads(args.result.read_text(encoding="utf-8"))
    selected = result.get("selected_ideal_upper_bound")
    if selected is None:
        selected = result.get("selected_screening_design", result.get("selected"))
    name = args.geometry or selected["geometry"]
    design = result.get("designs", {}).get(name)
    if design is None:
        geometry = result.get("geometries", {}).get(name, {})
        design = geometry.get("design")

    fig, axis = plt.subplots(figsize=(8.3, 8.0), constrained_layout=True)
    fig.patch.set_facecolor("white")
    for low, high, front, group in egg31_rings():
        color = "#F28E2B" if group == "original" else "#E15759"
        for x_sign in (-1.0, 1.0):
            axis.add_patch(Polygon(
                egg_polygon(low, high, front, x_sign), closed=True,
                facecolor=color, edgecolor="white", linewidth=0.55,
            ))
    if design is not None:
        draw_endcap(axis, design, "upstream")
        draw_endcap(axis, design, "downstream")

    axis.plot([], [], color="#F28E2B", linewidth=8, label="published 22 rings")
    axis.plot([], [], color="#E15759", linewidth=8,
              label="added 5 forward + 4 backward rings")
    if design is not None:
        axis.plot([], [], color="#7F7F7F", linewidth=5, label="PC lead")
        axis.plot([], [], color="#F2CF5B", linewidth=5, label="PC plastic")
    axis.scatter([0], [0], marker="x", color="black", s=36, zorder=20)
    axis.annotate("target", (0, 0), xytext=(5, 5), textcoords="offset points")
    axis.axhline(Z_OFFSET_CM, color="#666666", linewidth=0.7, linestyle=":")
    axis.text(43, Z_OFFSET_CM + 1, "BGOegg center z=-10 cm", fontsize=8,
              ha="right", color="#555555")
    axis.set_aspect("equal")
    axis.set_xlim(-65, 65)
    axis.set_ylim(-65, 65)
    axis.set_xlabel("transverse radius x (cm)")
    axis.set_ylabel("beam axis z (cm)")
    axis.grid(alpha=0.18)
    axis.legend(loc="upper right", fontsize=8, framealpha=0.95)
    pc_text = "no Photon Counter selected"
    if design is not None:
        pc_text = (
            f"PC={name}: {design['photon_counter']}, "
            f"{design['pc_n_layers']} x "
            f"({design['pc_pb_thickness_mm']:g} mm Pb + "
            f"{design['pc_scinti_thickness_mm']:g} mm plastic), "
            f"|z|={design['pc_z_front_cm']:g} cm"
        )
    axis.text(
        0.02, 0.02,
        "31 rings x 60 sectors; published frusta; BGO threshold = 3 MeV\n"
        + pc_text + "\n"
        "ideal active solids only: no support, PMT, services, pile-up, or real background",
        transform=axis.transAxes, fontsize=8.5, va="bottom",
        bbox={"facecolor": "white", "alpha": 0.88, "edgecolor": "none"},
    )
    fig.suptitle(
        "Exact extended BGOegg + Photon Counter diagnostic upper bound\n"
        "cross-section through crystal centers; target and z offset shown explicitly",
        fontsize=13,
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
