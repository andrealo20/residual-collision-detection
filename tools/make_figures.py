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

The figures are SVG with text converted to outlines, so they look the same on
every machine and need no font installed to view. The typeface is STIX, which
ships with matplotlib, so a regeneration on CI matches one on a laptop.

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
from matplotlib.text import Text  # noqa: E402
from matplotlib.ticker import FixedLocator, NullLocator  # noqa: E402

PAPER = "#ffffff"
INK = "#1d2330"
SOFT = "#5b6474"      # secondary text
FAINT = "#a3a9b4"     # reference marks: applied torque, guides
RULE = "#e6e8ec"      # hairlines
BAND = "#f1f2f4"      # contact window
BLUE = "#2a6496"      # the residual, and anything that is the method working
RED = "#c8502e"       # what model and sensor error do to it

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "docs", "data")
SWEEP_JSON = os.path.join(DATA, "sweep.json")

plt.rcParams.update({
    "font.family": "STIXGeneral",
    "mathtext.fontset": "stix",
    "font.size": 10,
    "text.color": INK,
    "axes.labelcolor": SOFT,
    "axes.labelsize": 9.5,
    "xtick.color": SOFT,
    "ytick.color": SOFT,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "axes.edgecolor": FAINT,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "svg.fonttype": "path",
    "svg.hashsalt": "qcollide",   # stable ids: unchanged figures diff clean
})


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


def find_case(sweep, prefix):
    for case in sweep["cases"]:
        if case["title"].lower().startswith(prefix):
            return case
    raise SystemExit(f"no case titled '{prefix}...' in the sweep")


def plain(ax):
    """Left and bottom rules only, no grid."""
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.tick_params(length=3, pad=3)


def dot_axis(ax):
    """Log x axis in Nm for the dot charts; rows carry their own labels."""
    ax.set_xscale("log")
    lo, hi = ax.get_xlim()
    ticks = [t for t in (0.001, 0.01, 0.1, 1, 10, 100) if lo <= t <= hi]
    ax.xaxis.set_major_locator(FixedLocator(ticks))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_xticklabels([f"{t:g}" for t in ticks])
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.yaxis.set_visible(False)
    for t in ticks:
        ax.axvline(t, color=RULE, linewidth=0.6, zorder=0)


def row_label(ax, y, text):
    """Row name in the left margin, clear of the data."""
    ax.text(-0.015, y, text, transform=ax.get_yaxis_transform(), color=INK,
            fontsize=9.5, ha="right", va="center")


def panel_title(ax, letter, text):
    ax.text(0, 1.06, letter, transform=ax.transAxes, fontsize=11,
            fontweight="bold", ha="left", va="bottom")
    ax.text(0.03, 1.06, text, transform=ax.transAxes, fontsize=10.5,
            ha="left", va="bottom")


def note(ax, x, y, text, **kw):
    """In-plot text on a paper patch, so gridlines pass behind it."""
    kw.setdefault("fontsize", 9)
    kw.setdefault("ha", "left")
    kw.setdefault("va", "center")
    return ax.text(x, y, text, zorder=5,
                   bbox=dict(facecolor=PAPER, edgecolor="none", pad=1.0), **kw)


def check_overlaps(fig, name):
    """Fail if any text touches other text or a plotted line or marker.

    Gridlines may pass behind text that sits on a paper patch.
    """
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    texts = [(t, t.get_window_extent(r).padded(1))
             for t in fig.findobj(Text)
             if t.get_visible() and t.get_text().strip()]
    problems = []
    for i, (a, ba) in enumerate(texts):
        for b, bb in texts[i + 1:]:
            if ba.overlaps(bb):
                problems.append(f"'{a.get_text()}' / '{b.get_text()}'")
    for ax in fig.axes:
        for ln in ax.lines:
            masked = ln.get_color() == RULE
            path = ln.get_transform().transform_path(ln.get_path())
            radius = r.points_to_pixels(ln.get_markersize()) / 2
            for t, box in texts:
                if masked and t.get_bbox_patch() is not None:
                    continue
                if ln.get_linestyle() not in ("None", ""):
                    hit = path.intersects_bbox(box, filled=False)
                else:
                    grown = box.padded(radius)
                    hit = any(grown.contains(x, y) for x, y in path.vertices)
                if hit:
                    problems.append(f"'{t.get_text()}' / a plotted line")
    if problems:
        raise SystemExit(f"{name}: overlapping elements:\n  "
                         + "\n  ".join(problems))


def save(fig, outdir, name):
    check_overlaps(fig, name)
    path = os.path.join(outdir, name)
    fig.savefig(path, bbox_inches="tight", pad_inches=0.12,
                facecolor=PAPER, metadata={"Date": None})
    plt.close(fig)
    return path


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


def fig_trace(sweep, outdir):
    """The residual following a contact, with and without model error."""
    gain = sweep["gain"]
    everything = lowest_usable(find_case(sweep, "everything"))
    if everything is None:
        raise SystemExit("the 'everything' case has no usable threshold")
    threshold = everything[0]

    fig, axes = plt.subplots(2, 1, figsize=(8.0, 5.2), sharex=True)
    fig.subplots_adjust(hspace=0.42)

    for ax, name, letter, title in (
        (axes[0], "trace_ideal.csv", "a", "Perfect model, ideal sensors"),
        (axes[1], "trace_realistic.csv", "b",
         "Friction, 5% parameter error, real sensors"),
    ):
        t, mag, ext = load_trace(name)
        inside = ext > 0
        t0, t1 = t[inside].min(), t[inside].max()
        peak = ext.max()

        ax.axvspan(t0, t1, color=BAND, linewidth=0, zorder=0)
        ax.plot(t, ext, color=FAINT, linewidth=1.1, zorder=2)
        ax.plot(t, mag, color=BLUE, linewidth=0.7, zorder=3)
        ax.text(t1 + 0.04, peak, f"applied, {peak:g} Nm", color=SOFT,
                fontsize=9, ha="left", va="center")

        ax.set_ylim(-0.4, 7.6)
        ax.set_yticks([0, 2, 4, 6])
        ax.set_ylabel("torque (Nm)")
        panel_title(ax, letter, title)
        plain(ax)

    # a: the step response is the observer's whole design, one time constant.
    ax = axes[0]
    tau_ms = 1000.0 / gain
    ax.text(2.05, 2.6, "residual, rises and decays with\n"
            f"time constant $1/K$ = {tau_ms:g} ms", color=BLUE, fontsize=9,
            ha="left", va="center")

    # b: the floor without contact, and the threshold that has to clear it.
    ax = axes[1]
    t, mag, ext = load_trace("trace_realistic.csv")
    # Before the contact only: the decay after it is the observer catching up.
    floor = mag[t < t[ext > 0].min()].max()
    ax.axhline(threshold, color=RED, linewidth=0.8, zorder=2)
    ax.text(0.08, threshold + 0.15,
            f"lowest threshold that stays silent, {threshold:g} Nm",
            color=RED, fontsize=9, ha="left", va="bottom")
    ax.text(0.08, floor + 0.3,
            f"residual with nothing touching the arm, up to {floor:.1f} Nm",
            color=BLUE, fontsize=9, ha="left", va="bottom")

    axes[1].set_xlabel("time (s)")
    axes[1].set_xlim(0, 3.0)
    return save(fig, outdir, "residual_trace.svg")


def fig_floors(sweep, outdir):
    """Where the residual floor comes from, every source on one scale."""
    baseline = sweep["friction_floors"][0]
    rows = [("Model error, ideal sensors", None, None)]
    rows.append(("none", baseline["floor"], False))
    for r in sweep["friction_floors"]:
        if r["scale"] > 0:
            rows.append((f"friction, {r['scale'] * 100:g}% of typical",
                         r["floor"], False))
    for r in sweep["parameter_floors"]:
        if r["error"] > 0:
            rows.append((f"masses {r['error'] * 100:g}% heavy",
                         r["floor"], False))
    rows.append(("Sensors, perfect model", None, None))
    sensors = sweep["sensor_floors"]
    worst = max(sensors, key=lambda r: r["floor"])
    for r in sensors:
        rows.append((r["case"], r["floor"], r is worst))

    fig, ax = plt.subplots(figsize=(5.6, 0.27 * len(rows) + 0.5))
    ax.set_xscale("log")
    ax.set_xlim(0.0006, 150)
    dot_axis(ax)
    ax.set_ylim(len(rows) - 0.4, -0.8)

    filtered = next(r for r in sensors if "hz filter" in r["case"].lower())
    for y, (label, value, highlight) in enumerate(rows):
        if value is None:
            ax.text(-0.015, y + 0.1, label, transform=ax.get_yaxis_transform(),
                    fontsize=10, fontweight="bold", ha="right", va="bottom")
            continue
        colour = RED if highlight else BLUE
        ax.plot([0.0006, value], [y, y], color=RULE, linewidth=0.6, zorder=1)
        ax.plot(value, y, "o", color=colour, markersize=5.5, zorder=3,
                markeredgecolor=PAPER, markeredgewidth=1.2)
        text = f"{value:.3f}"
        # The one comparison the sensor rows exist to make.
        if label == filtered["case"]:
            text += f", {worst['floor'] / value:.0f} times lower"
        note(ax, value * 1.4, y, text, color=SOFT, fontsize=8.5)
        row_label(ax, y, label)

    ax.set_xlabel("peak residual with no contact (Nm, log scale)")
    return save(fig, outdir, "residual_floor.svg")


def fig_tradeoff(sweep, outdir):
    """Lowest usable threshold and the smallest contact it still catches."""
    rows = []
    for case in sweep["cases"]:
        found = lowest_usable(case)
        if found is None:
            continue
        label = case["title"].replace(" (baseline)", "")
        rows.append((label, found[0], found[1]))
    if not rows:
        raise SystemExit("no case in the sweep has a usable threshold")

    fig, ax = plt.subplots(figsize=(4.6, 0.36 * len(rows) + 0.9))
    ax.set_xscale("log")
    ax.set_xlim(0.05, 20)
    dot_axis(ax)
    best = rows[0][2]
    ax.set_ylim(len(rows) - 0.4, -0.6)
    for y, (label, thr, mdt) in enumerate(rows):
        ax.plot([0.05, min(thr, mdt)], [y, y], color=RULE, linewidth=0.6,
                zorder=1)
        ax.plot([thr, mdt], [y, y], color=SOFT, linewidth=0.8, zorder=2)
        ax.plot(thr, y, "o", markersize=6, markerfacecolor=PAPER,
                markeredgecolor=SOFT, markeredgewidth=1.0, zorder=3)
        ax.plot(mdt, y, "o", markersize=6, color=RED if y else BLUE,
                markeredgecolor=PAPER, markeredgewidth=1.2, zorder=4)
        text = f"{mdt:.2f} Nm"
        if y == len(rows) - 1:
            text += f", {mdt / best:.0f}× the baseline"
        note(ax, max(thr, mdt) * 1.3, y, text, color=INK)
        row_label(ax, y, label)

    # Key, written out above the plot rather than boxed inside it.
    key = ax.transAxes
    ax.plot(0.01, 1.13, "o", markersize=6, markerfacecolor=PAPER,
            markeredgecolor=SOFT, markeredgewidth=1.0, clip_on=False,
            transform=key)
    ax.text(0.04, 1.13, "lowest threshold that stays silent",
            color=SOFT, fontsize=9, va="center", transform=key)
    ax.plot(0.01, 1.05, "o", markersize=6, color=RED,
            markeredgecolor=PAPER, markeredgewidth=1.2, clip_on=False,
            transform=key)
    ax.text(0.04, 1.05, "smallest contact it still catches",
            color=SOFT, fontsize=9, va="center", transform=key)

    ax.set_xlabel("torque (Nm, log scale)")
    return save(fig, outdir, "tradeoff.svg")


def main():
    outdir = (sys.argv[1] if len(sys.argv) > 1
              else os.path.join(ROOT, "docs", "images"))
    os.makedirs(outdir, exist_ok=True)
    sweep = load_sweep()

    print("wrote", fig_trace(sweep, outdir))
    print("wrote", fig_floors(sweep, outdir))
    print("wrote", fig_tradeoff(sweep, outdir))
    return 0


if __name__ == "__main__":
    sys.exit(main())
