#!/usr/bin/env python3
"""
Check the README's measurement tables against the sweep that produced them.

The README quotes residual floors and threshold trade-offs. Nothing used to
connect those tables to the program that measures them, so they could drift
apart silently, and the figures alongside them were drawn from numbers typed
into make_figures.py by hand.

The rule enforced here is deliberately simple, because it survives the tables
being reorganised: every measured quantity the README states must appear in
the sweep output. A number that appears in the prose but in no measurement is
either stale or invented, and both are worth failing over.

Usage:
    ./build/bench/qcollide_sweep --json docs/data/sweep.json
    python3 tools/check_readme.py [sweep.json] [README.md]
"""

import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

# The README rounds floors to three decimals and torques to two, so a quoted
# value is accepted when it matches a measurement to within half a unit of the
# last digit it shows.
def close(quoted, measured, decimals):
    return abs(quoted - measured) <= 0.5 * 10 ** (-decimals) + 1e-12


def collect(sweep):
    """Every number the sweep measured, grouped by what it is."""
    floors = [r["floor"] for r in sweep["friction_floors"]]
    floors += [r["floor"] for r in sweep["parameter_floors"]]
    sensor = [r["floor"] for r in sweep["sensor_floors"]]

    thresholds, detectables = set(), []
    for case in sweep["cases"]:
        for row in case["thresholds"]:
            thresholds.add(round(row["threshold"], 4))
            if row["usable"] and row["min_detectable"] is not None:
                detectables.append(row["min_detectable"])
    return floors, sensor, sorted(thresholds), detectables


def table_rows(text, header_words):
    """Rows of the first markdown table whose header holds every given word."""
    rows, in_table, matched = [], False, False
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped.startswith("|"):
            if matched and in_table:
                break
            in_table = False
            continue
        if set(stripped) <= set("|-: "):
            continue
        cells = [c.strip().strip("*` ") for c in stripped.strip("|").split("|")]
        if not in_table:
            in_table = True
            matched = all(w.lower() in stripped.lower() for w in header_words)
            if matched:
                continue
        if matched:
            rows.append(cells)
    return rows


def numbers_in(cell):
    return [float(x) for x in re.findall(r"-?\d+\.\d+|-?\d+", cell)]


def main():
    sweep_path = (sys.argv[1] if len(sys.argv) > 1
                  else os.path.join(ROOT, "docs", "data", "sweep.json"))
    readme_path = (sys.argv[2] if len(sys.argv) > 2
                   else os.path.join(ROOT, "README.md"))

    if not os.path.exists(sweep_path):
        raise SystemExit(
            f"{sweep_path} not found. Produce it first:\n"
            "    ./build/bench/qcollide_sweep --json docs/data/sweep.json")

    with open(sweep_path, encoding="utf-8") as fh:
        sweep = json.load(fh)
    with open(readme_path, encoding="utf-8") as fh:
        text = fh.read()

    floors, sensor_floors, thresholds, detectables = collect(sweep)
    all_floors = floors + sensor_floors
    problems = []

    # Residual floor tables: one number per row, in Nm.
    for header in (["model error", "peak residual"],
                   ["measurement condition", "peak residual"],
                   ["measurement", "peak residual"]):
        for cells in table_rows(text, header):
            if len(cells) < 2:
                continue
            quoted = numbers_in(cells[1])
            if not quoted:
                continue
            value = quoted[0]
            if not any(close(value, m, 3) for m in all_floors):
                problems.append(
                    f"residual floor {value} Nm (row '{cells[0]}') "
                    f"appears in no measurement")

    # Trade-off tables: a threshold and the smallest contact detected there.
    for header in (["lowest usable threshold", "smallest contact"],
                   ["lowest usable threshold"]):
        for cells in table_rows(text, header):
            if len(cells) < 3:
                continue
            thr = numbers_in(cells[1])
            mdt = numbers_in(cells[2])
            if thr and not any(close(thr[0], t, 2) for t in thresholds):
                problems.append(
                    f"threshold {thr[0]} Nm (row '{cells[0]}') is not one the "
                    f"sweep tries: {thresholds}")
            if mdt and not any(close(mdt[0], d, 2) for d in detectables):
                problems.append(
                    f"minimum detectable torque {mdt[0]} Nm (row '{cells[0]}') "
                    f"appears in no measurement")

    # Several header patterns can match the same table, so the same row can be
    # reported twice. Order is kept because it follows the README.
    seen = set()
    problems = [p for p in problems if not (p in seen or seen.add(p))]

    if problems:
        print(f"README.md disagrees with {os.path.basename(sweep_path)}: "
              f"{len(problems)} problems")
        for p in problems:
            print(f"  {p}")
        return 1

    print(f"README.md agrees with {os.path.basename(sweep_path)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
