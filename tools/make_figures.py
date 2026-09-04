#!/usr/bin/env python3
"""
Regenerate the figures used in the README.

Every number plotted is read from a measurement: the traces from
bench/qcollide_trace, the residual floors and the threshold trade-off from the
JSON that bench/qcollide_sweep writes with --json.

That last part used to be a lie. The floor and trade-off figures were drawn
from lists of numbers transcribed into this file by hand, so re-running the
sweep changed the tables and left the figures showing whatever had been typed
in last, while the README claimed the figures could not drift from the numbers
they illustrate. There are no constants here now, and a missing input is an
error rather than a fallback.

Usage:
    ./build/bench/qcollide_sweep --json docs/data/sweep.json
    ./build/bench/qcollide_trace 0.0 0.00 0 > docs/data/trace_ideal.csv
    ./build/bench/qcollide_trace 1.0 0.05 1 > docs/data/trace_realistic.csv
    python3 tools/make_figures.py [output_dir]

Requires: numpy, matplotlib.
"""

import csv
import json
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
SWEEP_JSON = os.path.join(DATA, "sweep.json")


def load_sweep():
    """Read the sweep results, or say how to produce them."""
    if not os.path.exists(SWEEP_JSON):
        raise SystemExit(
            f"{SWEEP_JSON} not found. Produce it first:\n"
            "    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release\n"
            "    cmake --build build --parallel\n"
            "    ./build/bench/qcollide_sweep --json docs/data/sweep.json\n"
        )
    with open(SWEEP_JSON, encoding="utf-8") as fh:
        return json.load(fh)


def lowest_usable(case):
    """The lowest threshold that does not fire on model error alone.

    Returns (threshold, smallest contact detected there), or None when every
    threshold in the sweep either chatters or misses everything.
    """
    for row in sorted(case["thresholds"], key=lambda r: r["threshold"]):
        if row["usable"] and row["min_detectable"] is not None:
            return row["threshold"], row["min_detectable"]
    return None


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
    with open(path, encoding="utf-8") as f:
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
        (axes[0], "trace_ideal.csv", "Perfect model, ideal sensors"),
        (axes[1], "trace_realistic.csv",
         "Friction, 5% parameter error, real sensors"),
    ):
        t, mag, ext = load_trace(name)

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


def fig_floors(sweep, outdir):
    """Where the residual floor comes from."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.0))
    fig.subplots_adjust(wspace=0.32)

    fx = [r["scale"] * 100 for r in sweep["friction_floors"]]
    fy = [r["floor"] for r in sweep["friction_floors"]]
    px = [r["error"] * 100 for r in sweep["parameter_floors"]]
    py = [r["floor"] for r in sweep["parameter_floors"]]

    ax1.plot(fx, fy, "o-", color=ACCENT, linewidth=2, markersize=6,
             zorder=3, label="friction (% of typical)")
    ax1.plot(px, py, "s-", color=DARK, linewidth=2, markersize=5,
             zorder=3, label="inertia error (%)")
    ax1.set_xlabel("severity (%)")
    ax1.set_xlim(-4, max(fx + px) * 1.04)
    ax1.set_ylabel("residual floor (Nm)")
    ax1.set_title("Model error", color=INK, fontsize=11, loc="left", pad=8)
    ax1.legend(frameon=False, fontsize=9)
    style(ax1)

    # Log scale: the raw difference is two orders of magnitude above
    # everything else, and a linear axis would hide the rest entirely.
    rows = sweep["sensor_floors"]
    labels = [r["case"].replace(", ", ",\n") for r in rows]
    vals = [r["floor"] for r in rows]
    worst = max(range(len(vals)), key=lambda i: vals[i])
    colours = [ACCENT if i == worst else DARK for i in range(len(vals))]

    ax2.bar(range(len(vals)), vals, color=colours, zorder=3, width=0.6)
    ax2.set_yscale("log")
    ax2.set_xticks(range(len(vals)))
    ax2.set_xticklabels(labels, fontsize=7.5)
    ax2.set_ylabel("residual floor (Nm), log scale")
    ax2.set_title("Sensors", color=INK, fontsize=11, loc="left", pad=8)
    for i, v in enumerate(vals):
        ax2.annotate(f"{v:.3f}", (i, v), textcoords="offset points",
                     xytext=(0, 4), ha="center", fontsize=8, color=INK,
                     zorder=4)
    ax2.set_ylim(min(v for v in vals if v > 0) * 0.3, max(vals) * 8)
    style(ax2)

    fig.suptitle("What raises the residual when nothing is touching the arm",
                 color=MUTED, fontsize=9, y=1.02, x=0.005, ha="left")
    fig.tight_layout(w_pad=2.5)
    path = os.path.join(outdir, "residual_floor.png")
    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return path


def fig_tradeoff(sweep, outdir):
    """Lowest usable threshold and the smallest contact it still catches."""
    names, thr, mdt = [], [], []
    for case in sweep["cases"]:
        found = lowest_usable(case)
        if found is None:
            continue
        names.append(case["title"].replace(", ", ",\n").replace(": ", ":\n"))
        thr.append(found[0])
        mdt.append(found[1])

    if not names:
        raise SystemExit("no case in the sweep has a usable threshold")

    fig, ax = plt.subplots(figsize=(8.6, 4.4))
    x = np.arange(len(names))
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
    ax.set_xticklabels(names, fontsize=7.5)
    ax.set_ylabel("torque (Nm)")
    ax.set_ylim(0, max(thr + mdt) * 1.25)
    ax.set_title("Detection sensitivity against realism", color=INK,
                 fontsize=11, loc="left", pad=8)
    ax.legend(frameon=False, fontsize=9, loc="upper left")
    style(ax)

    best = min(mdt)
    worst = max(mdt)
    fig.suptitle(f"A {best:.2f} Nm touch becomes a {worst:.2f} Nm push once "
                 "ordinary model error and sensors are present",
                 color=MUTED, fontsize=9, y=1.02, x=0.005, ha="left")
    fig.tight_layout()
    path = os.path.join(outdir, "tradeoff.png")
    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return path


def main():
    outdir = (sys.argv[1] if len(sys.argv) > 1
              else os.path.join(ROOT, "docs", "images"))
    os.makedirs(outdir, exist_ok=True)
    sweep = load_sweep()

    print("wrote", fig_trace(outdir))
    print("wrote", fig_floors(sweep, outdir))
    print("wrote", fig_tradeoff(sweep, outdir))
    return 0


if __name__ == "__main__":
    sys.exit(main())
