"""Prepare point-cloud samples for the finite-difference reconstruction code.

The solver expects:

* nt_points.m: one integer, the number of point samples.
* fun_data.m: x, y, z coordinates written one value per line.

This utility reads a plain text point-cloud file with either three columns
(`x y z`) or solver-style one-value-per-line triples, optionally normalizes
the samples into a unit box, and writes the solver input files.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


def read_points(path: Path) -> np.ndarray:
    values = np.loadtxt(path, dtype=float)
    if values.ndim == 1:
        if values.size % 3 != 0:
            raise ValueError("one-column point data must contain x/y/z triples")
        points = values.reshape((-1, 3))
    elif values.ndim == 2 and values.shape[1] >= 3:
        points = values[:, :3]
    else:
        raise ValueError("point data must have three columns or one-value-per-line triples")

    if points.size == 0:
        raise ValueError("point data is empty")
    return points


def normalize_points(points: np.ndarray, scaling_factor: float) -> np.ndarray:
    lower = points.min(axis=0)
    upper = points.max(axis=0)
    span = np.max(upper - lower)
    if span <= 0:
        raise ValueError("point data has zero spatial extent")

    center = (upper + lower) / 2.0
    return (points - center) * (scaling_factor / span) + np.array([0.5, 0.5, 0.5])


def write_solver_inputs(points: np.ndarray, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    count_file = output_dir / "nt_points.m"
    data_file = output_dir / "fun_data.m"

    count_file.write_text(f"{len(points)}\n", encoding="utf-8")
    with data_file.open("w", encoding="utf-8") as f:
        for x, y, z in points:
            f.write(f"{x:.10f}\n{y:.10f}\n{z:.10f}\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path, help="plain text point-cloud input")
    parser.add_argument("--output-dir", required=True, type=Path, help="folder for solver inputs")
    parser.add_argument("--normalize", action="store_true", help="normalize points into a unit box")
    parser.add_argument("--scale", default=0.8, type=float, help="scale factor used with --normalize")
    args = parser.parse_args()

    points = read_points(args.input)
    if args.normalize:
        points = normalize_points(points, args.scale)
    write_solver_inputs(points, args.output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
