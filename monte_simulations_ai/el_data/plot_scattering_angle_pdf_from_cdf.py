#!/usr/bin/env python3
"""Build and plot scattering-angle PDFs from inverse-CDF tables stored in a CDL/CDF file."""

from __future__ import annotations

import argparse
import os
import tempfile
from pathlib import Path

THIS_DIR = Path(__file__).resolve().parent
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "matplotlib"))

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

from cdf_compat import ElasticCdf, probability_axis, reshape_scattering_angles


DEFAULT_CDF = THIS_DIR / "Krstic" / "krstic_dd_total_elastic_integral_priority.cdf"
DEFAULT_OUT_DIR = THIS_DIR / "figure" / "scattering_angle_pdf"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cdf",
        type=Path,
        default=DEFAULT_CDF,
        help="Input CDL/CDF file that contains scattering_angle data.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=DEFAULT_OUT_DIR,
        help="Directory for figure outputs.",
    )
    parser.add_argument(
        "--sample-energies",
        type=float,
        nargs="*",
        default=None,
        help="Requested sample energies [eV/amu] for line plots. Nearest grid points are used.",
    )
    parser.add_argument(
        "--num-samples",
        type=int,
        default=6,
        help="Number of representative energies to use when --sample-energies is omitted.",
    )
    parser.add_argument(
        "--theta-grid-points",
        type=int,
        default=512,
        help="Number of theta points for the interpolated PDF heatmap.",
    )
    return parser.parse_args()


def build_pdf_from_inverse_cdf(
    prob_axis: np.ndarray,
    theta_values: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, float]:
    prob = np.asarray(prob_axis, dtype=float)
    theta = np.asarray(theta_values, dtype=float)
    if prob.ndim != 1 or theta.ndim != 1:
        raise ValueError("Probability axis and theta values must both be 1D arrays.")
    if prob.size != theta.size:
        raise ValueError("Probability axis and theta values must have the same length.")
    if prob.size < 2:
        raise ValueError("At least two inverse-CDF samples are required.")

    delta_prob = np.diff(prob)
    delta_theta = np.diff(theta)
    width = np.abs(delta_theta)
    if np.any(width <= 0.0):
        raise ValueError("Theta table must be strictly monotonic to build a PDF.")

    theta_mid = 0.5 * (theta[:-1] + theta[1:])
    pdf = np.abs(delta_prob / delta_theta)
    normalization = float(np.sum(pdf * width))

    if theta_mid[0] > theta_mid[-1]:
        return theta_mid[::-1], pdf[::-1], normalization
    return theta_mid, pdf, normalization


def choose_sample_indices(
    energy_axis: np.ndarray,
    requested_energies: list[float] | None,
    num_samples: int,
) -> list[int]:
    energy = np.asarray(energy_axis, dtype=float)
    if requested_energies:
        selected: list[int] = []
        for requested in requested_energies:
            index = int(np.argmin(np.abs(energy - float(requested))))
            if index not in selected:
                selected.append(index)
        return selected

    count = max(1, min(int(num_samples), energy.size))
    return list(np.linspace(0, energy.size - 1, count, dtype=int))


def plot_sample_lines(
    theta_rows: list[np.ndarray],
    pdf_rows: list[np.ndarray],
    energy_axis: np.ndarray,
    sample_indices: list[int],
    output_path: Path,
) -> None:
    fig, ax = plt.subplots(figsize=(9.0, 5.8))
    for index in sample_indices:
        ax.plot(
            theta_rows[index],
            pdf_rows[index],
            lw=1.8,
            label=f"E = {energy_axis[index]:.4g} eV/amu",
        )

    ax.set_xlim(0.0, float(np.pi))
    ax.set_yscale("log")
    ax.set_xlabel("scattering angle θ [rad]")
    ax.set_ylabel("PDF p(θ) [rad⁻¹]")
    ax.set_title("Scattering-angle PDF reconstructed from inverse CDF")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def plot_heatmap(
    theta_rows: list[np.ndarray],
    pdf_rows: list[np.ndarray],
    energy_axis: np.ndarray,
    theta_grid_points: int,
    output_path: Path,
) -> None:
    theta_common = np.linspace(0.0, float(np.pi), int(theta_grid_points))
    pdf_grid = np.zeros((len(pdf_rows), theta_common.size), dtype=float)
    for row_index, (theta_row, pdf_row) in enumerate(zip(theta_rows, pdf_rows, strict=True)):
        pdf_grid[row_index] = np.interp(theta_common, theta_row, pdf_row, left=0.0, right=0.0)

    positive = pdf_grid[pdf_grid > 0.0]
    if positive.size == 0:
        raise ValueError("Interpolated PDF grid has no positive values.")

    fig, ax = plt.subplots(figsize=(9.0, 5.8))
    mesh = ax.pcolormesh(
        theta_common,
        energy_axis,
        pdf_grid,
        shading="auto",
        cmap="viridis",
        norm=LogNorm(vmin=float(positive.min()), vmax=float(positive.max())),
    )
    ax.set_xlim(0.0, float(np.pi))
    ax.set_yscale("log")
    ax.set_xlabel("scattering angle θ [rad]")
    ax.set_ylabel("E [eV/amu]")
    ax.set_title("Scattering-angle PDF map reconstructed from inverse CDF")
    fig.colorbar(mesh, ax=ax, label="PDF p(θ) [rad⁻¹]")
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def main() -> None:
    args = parse_args()

    cdf = ElasticCdf.from_file(args.cdf)
    prob_axis = probability_axis(cdf)
    energy_axis = cdf.get_axis(5, 1)
    angle_rows = reshape_scattering_angles(cdf.get_block(5), cdf)

    theta_rows: list[np.ndarray] = []
    pdf_rows: list[np.ndarray] = []
    normalizations: list[float] = []
    for theta_values in angle_rows:
        theta_mid, pdf, normalization = build_pdf_from_inverse_cdf(prob_axis, theta_values)
        theta_rows.append(theta_mid)
        pdf_rows.append(pdf)
        normalizations.append(normalization)

    sample_indices = choose_sample_indices(energy_axis, args.sample_energies, args.num_samples)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    line_path = args.out_dir / "scattering_angle_pdf_samples.png"
    heatmap_path = args.out_dir / "scattering_angle_pdf_map.png"
    plot_sample_lines(theta_rows, pdf_rows, energy_axis, sample_indices, line_path)
    plot_heatmap(theta_rows, pdf_rows, energy_axis, args.theta_grid_points, heatmap_path)

    normalization_error = np.max(np.abs(np.asarray(normalizations) - 1.0))
    print(f"input cdf        : {args.cdf}")
    print(f"sample-line plot : {line_path}")
    print(f"heatmap plot     : {heatmap_path}")
    print(f"max |∫p(θ)dθ-1|   : {normalization_error:.3e}")


if __name__ == "__main__":
    main()
