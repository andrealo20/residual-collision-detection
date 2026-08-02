#!/usr/bin/env python3
"""
Regenerate the figures used in the README.

Every number plotted comes from a measurement rather than a model: the traces
from bench/qcollide_trace, the residual floors and threshold trade-offs from
bench/qcollide_sweep and the test suite. The script exists so the figures can
be rebuilt when the measurements change, and so that where each one came from
is written down rather than remembered.

Usage:
    ./build/bench/qcollide_trace 0.0 0.00 0 > docs/data/trace_ideal.csv
    ./build/bench/qcollide_trace 1.0 0.05 1 > docs/data/trace_realistic.csv
    python3 tools/make_figures.py [output_dir]

Requires: numpy, matplotlib.
"""

import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402
import numpy as np                # noqa: E402

INK = "#1b1f24"
ACCENT = "#c2410c"
MUTED = "#94a3b8"
GRID = "#e2e8f0"
DARK = "#334155"

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "docs", "data")

# --- measured: bench/qcollide_sweep and tests/test_sensors ----------------
FLOOR_FRICTION = [(0.00, 0.001), (0.25, 0.603), (0.50, 1.101), (1.00, 1.866)]
FLOOR_PARAM = [(0.00, 0.001), (0.02, 0.178), (0.05, 0.445), (0.10, 0.890)]
FLOOR_SENSOR = [
    ("exact velocity", 0.003),
    ("+ torque noise", 0.010),
    ("raw encoder\ndifference", 26.313),
    ("+ 50 Hz filter", 0.477),
    ("+ torque noise", 0.480),
]

# --- measured: bench/qcollide_sweep --------------------------------------
CASES = [
    ("perfect model\nideal sensors", 0.10, 0.10),
    ("real sensors\nperfect model", 1.00, 0.88),
    ("friction only\nideal sensors", 2.00, 3.61),
    ("friction + 5%\nparameters", 4.00, 5.28),
    ("everything", 4.00, 4.80),
]


def style(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
    ax.tick_params(colors=INK, labelsize=9)
    ax.grid(True, color=GRID, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)


def load_trace(name):
    path = os.path.join(DATA, name)
    if not os.path.exists(path):
        raise SystemExit(f"{path} not found: run bench/qcollide_trace first")
    t, mag, ext = [], [], []
    with open(path) as f:
        for row in csv.DictReader(f):
            t.append(float(row["t"]))
            mag.append(float(row["magnitude"]))
            ext.append(float(row["tau_ext0"]))
    return np.array(t), np.array(mag), np.array(ext)


def fig_trace(outdir):
    """The residual following a contact, with and without model error."""
    fig, axes = plt.subplots(2, 1, figsize=(9.5, 5.8), sharex=True)
    # Extra vertical room: without it the lowest tick of the upper panel and
    # the highest tick of the lower one collide, which an automated check
    # caught and the eye did not.
    fig.subplots_adjust(hspace=0.35)

    for ax, name, title in (
        (axes[0], "trace_ideal.csv",
         "Perfect model, ideal sensors"),
        (axes[1], "trace_realistic.csv",
         "Friction, 5% parameter error, real sensors"),
    ):
        t, mag, ext = load_trace(name)

        # Shade the contact window, so the reader can see the residual rise
        # inside it rather than having to infer the timing.
        inside = ext > 0
        if inside.any():
            ax.axvspan(t[inside].min(), t[inside].max(), color=ACCENT,
                       alpha=0.10, zorder=1, label="contact")

        ax.plot(t, ext, color=MUTED, linewidth=1.4, linestyle="--",
                zorder=2, label="applied external torque")
        ax.plot(t, mag, color=ACCENT, linewidth=1.4, zorder=3,
                label="residual magnitude")
        ax.set_ylabel("torque (Nm)")
        ax.set_title(title, color=INK, fontsize=11, loc="left", pad=8)
        ax.set_ylim(-0.8, 7.8)
        style(ax)

    axes[0].legend(frameon=False, fontsize=9, loc="upper left")
    axes[1].set_xlabel("time (s)")

    fig.suptitle("The residual tracks the external torque, and model error "
                 "raises the floor it has to rise above",
                 color=MUTED, fontsize=9, y=1.01, x=0.005, ha="left")
    fig.tight_layout(h_pad=1.8)
    path = os.path.join(outdir, "residual_trace.png")
    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return path


def fig_floors(outdir):
    """Where the residual floor comes from."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.0))
    # Horizontal room for the same reason: the x-tick labels of the left
    # panel otherwise reach the y-tick labels of the right one.
    fig.subplots_adjust(wspace=0.32)

    # Model error, linear scale: the two sources are comparable in size.
    fx = [v * 100 for v, _ in FLOOR_FRICTION]
    fy = [f for _, f in FLOOR_FRICTION]
    px = [v * 100 for v, _ in FLOOR_PARAM]
    py = [f for _, f in FLOOR_PARAM]

    ax1.plot(fx, fy, "o-", color=ACCENT, linewidth=2, markersize=6,
             zorder=3, label="friction (% of typical)")
    ax1.plot(px, py, "s-", color=DARK, linewidth=2, markersize=5,
             zorder=3, label="inertia error (%)")
    ax1.set_xlabel("severity (%)")
    # Capped at the data: the automatic limit ran to 120 % and its tick label
    # reached the log-scale labels of the right-hand panel.
    ax1.set_xlim(-4, 104)
    ax1.set_ylabel("residual floor (Nm)")
    ax1.set_title("Model error", color=INK, fontsize=11, loc="left", pad=8)
    ax1.legend(frameon=False, fontsize=9)
    style(ax1)

    # Sensors, log scale: the raw difference is two orders of magnitude
    # above everything else, and a linear axis would hide the rest entirely.
    labels = [n for n, _ in FLOOR_SENSOR]
    vals = [v for _, v in FLOOR_SENSOR]
    colours = [DARK, DARK, ACCENT, DARK, DARK]
    ax2.bar(range(len(vals)), vals, color=colours, zorder=3, width=0.6)
    ax2.set_yscale("log")
    ax2.set_xticks(range(len(vals)))
    ax2.set_xticklabels(labels, fontsize=8)
    ax2.set_ylabel("residual floor (Nm), log scale")
    ax2.set_title("Sensors", color=INK, fontsize=11, loc="left", pad=8)
    for i, v in enumerate(vals):
        ax2.annotate(f"{v:.3f}", (i, v), textcoords="offset points",
                     xytext=(0, 4), ha="center", fontsize=8, color=INK,
                     zorder=4)
    ax2.set_ylim(1e-3, 200)
    style(ax2)

    fig.suptitle("What raises the residual when nothing is touching the arm",
                 color=MUTED, fontsize=9, y=1.02, x=0.005, ha="left")
    fig.tight_layout(w_pad=2.5)
    path = os.path.join(outdir, "residual_floor.png")
    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return path


def fig_tradeoff(outdir):
    """Lowest usable threshold and the smallest contact it still catches."""
    fig, ax = plt.subplots(figsize=(8.2, 4.2))

    names = [c[0] for c in CASES]
    thr = [c[1] for c in CASES]
    mdt = [c[2] for c in CASES]
    x = np.arange(len(CASES))
    width = 0.38

    b1 = ax.bar(x - width / 2, thr, width, color=DARK, zorder=3,
                label="lowest usable threshold")
    b2 = ax.bar(x + width / 2, mdt, width, color=ACCENT, zorder=3,
                label="smallest contact detected")

    for bars in (b1, b2):
        for r in bars:
            ax.annotate(f"{r.get_height():.2f}",
                        (r.get_x() + r.get_width() / 2, r.get_height()),
                        textcoords="offset points", xytext=(0, 3),
                        ha="center", fontsize=8, color=INK, zorder=4)

    ax.set_xticks(x)
    ax.set_xticklabels(names, fontsize=8.5)
    ax.set_ylabel("torque (Nm)")
    ax.set_ylim(0, 6.6)
    ax.set_title("Detection sensitivity against realism", color=INK,
                 fontsize=11, loc="left", pad=8)
    ax.legend(frameon=False, fontsize=9, loc="upper left")
    style(ax)

    fig.suptitle("A 0.1 Nm touch becomes a 4.8 Nm push once ordinary model "
                 "error and sensors are present",
                 color=MUTED, fontsize=9, y=1.02, x=0.005, ha="left")
    fig.tight_layout()
    path = os.path.join(outdir, "tradeoff.png")
    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return path


def main():
    outdir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "docs",
                                                                "images")
    os.makedirs(outdir, exist_ok=True)
    for fn in (fig_trace, fig_floors, fig_tradeoff):
        print("wrote", fn(outdir))
    return 0


if __name__ == "__main__":
    sys.exit(main())
