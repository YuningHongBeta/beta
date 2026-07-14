#!/usr/bin/env python3
"""Draw the simulated current BGO ball and exact extended BGOegg cross-sections."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Polygon


PHI_DEG = 6.0
FORWARD_B = 1.3
CRYSTAL_LENGTH_CM = 22.0
BALL_THETA_MIN_DEG = 5.666
BALL_THETA_MAX_DEG = 170.302


def forward_low_deg(crystal_type: int) -> float:
    delta = math.radians(PHI_DEG)
    return math.degrees(
        math.pi / 2
        - math.atan(FORWARD_B * math.tan((crystal_type + 1) * delta / FORWARD_B))
    )


def egg31_rings() -> list[tuple[float, float, float, str]]:
    rings = []
    for crystal_type in range(17, -1, -1):
        low = forward_low_deg(crystal_type)
        high = 90.0 if crystal_type == 0 else forward_low_deg(crystal_type - 1)
        middle = math.radians(0.5 * (low + high))
        radius = 20.0 / math.sqrt(
            math.sin(middle) ** 2 + 0.25 * math.cos(middle) ** 2
        )
        rings.append((low, high, radius, "extended" if crystal_type >= 13 else "original"))
    for crystal_type in range(13):
        low = 90.0 + PHI_DEG * crystal_type
        rings.append((low, low + PHI_DEG, 20.0,
                      "extended" if crystal_type >= 9 else "original"))
    return rings


def radial_polygon(theta_low: float, theta_high: float,
                   r_front: float, r_rear: float, sign: float) -> list[tuple[float, float]]:
    points = []
    for radius, theta in (
        (r_front, theta_low),
        (r_rear, theta_low),
        (r_rear, theta_high),
        (r_front, theta_high),
    ):
        angle = math.radians(theta)
        points.append((sign * radius * math.sin(angle), radius * math.cos(angle)))
    return points


def draw_ball(axis: plt.Axes) -> None:
    cos_min = math.cos(math.radians(BALL_THETA_MIN_DEG))
    cos_max = math.cos(math.radians(BALL_THETA_MAX_DEG))
    boundaries = [
        math.degrees(math.acos(cos_min + (cos_max - cos_min) * index / 20.0))
        for index in range(21)
    ]
    for low, high in zip(boundaries[:-1], boundaries[1:]):
        for sign in (-1.0, 1.0):
            axis.add_patch(Polygon(
                radial_polygon(low, high, 30.0, 50.0, sign),
                closed=True, facecolor="#4C78A8", edgecolor="white", linewidth=0.7,
            ))
    axis.text(0.02, 0.02,
              "20 rings x 40 sectors = 800 cells\n"
              "equal-solid-angle segmentation\n"
              "R = 30--50 cm; polar range = 5.666--170.302 deg",
              transform=axis.transAxes, fontsize=9, va="bottom")
    axis.set_title("Current BGO ball (a20x40)")


def draw_egg(axis: plt.Axes) -> None:
    for low, high, front, group in egg31_rings():
        color = "#F28E2B" if group == "original" else "#E15759"
        for sign in (-1.0, 1.0):
            axis.add_patch(Polygon(
                radial_polygon(low, high, front, front + CRYSTAL_LENGTH_CM, sign),
                closed=True, facecolor=color, edgecolor="white", linewidth=0.7,
            ))
    axis.plot([], [], color="#F28E2B", linewidth=8, label="published 22 rings")
    axis.plot([], [], color="#E15759", linewidth=8,
              label="added 5 forward + 4 backward rings")
    axis.legend(loc="upper right", fontsize=8, frameon=True,
                facecolor="white", framealpha=0.92)
    axis.text(0.02, 0.02,
              "31 rings x 60 sectors = 1860 cells\n"
              "published frustum rule; crystal length = 22 cm\n"
              "polar range = 5.336--168 deg; no gaps/support/PC",
              transform=axis.transAxes, fontsize=9, va="bottom")
    axis.set_title("Extended BGOegg (exact frusta)")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdf", type=Path, required=True)
    parser.add_argument("--png", type=Path, required=True)
    args = parser.parse_args()

    fig, axes = plt.subplots(1, 2, figsize=(12.0, 6.2), constrained_layout=True)
    fig.patch.set_facecolor("white")
    for axis, draw in zip(axes, (draw_ball, draw_egg)):
        draw(axis)
        axis.scatter([0], [0], marker="x", color="black", s=30, zorder=10)
        axis.annotate("target", (0, 0), xytext=(4, 5), textcoords="offset points",
                      fontsize=8)
        axis.set_aspect("equal")
        axis.set_xlim(-65, 65)
        axis.set_ylim(-65, 65)
        axis.set_xlabel("transverse radius x (cm)")
        axis.set_ylabel("beam axis z (cm)")
        axis.grid(alpha=0.18)
    fig.suptitle(
        "Geometry punch drawing: current BGO ball vs segment-preserving extended BGOegg\n"
        "Geant4 model used for the 2026-07-14 3 MeV study; cross-section through cell centers",
        fontsize=13,
    )
    args.pdf.parent.mkdir(parents=True, exist_ok=True)
    args.png.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(args.pdf) as pdf:
        pdf.savefig(fig, bbox_inches="tight", facecolor="white")
    fig.savefig(args.png, dpi=180, bbox_inches="tight",
                facecolor="white", transparent=False)
    plt.close(fig)
    print(args.pdf)
    print(args.png)


if __name__ == "__main__":
    main()
