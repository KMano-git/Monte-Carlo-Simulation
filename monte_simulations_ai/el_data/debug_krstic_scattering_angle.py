#!/usr/bin/env python3
"""Debug Krstic scattering-angle tables against the repository Bachmann baseline."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def load_path_module(module_name: str, module_path: Path):
    spec = importlib.util.spec_from_file_location(module_name, module_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Failed to load module {module_name!r} from {module_path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def parse_args() -> argparse.Namespace:
    this_dir = Path(__file__).resolve().parent
    krstic_dir = this_dir / "Krstic"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--bachmann-cdf",
        default=str(this_dir / "Bachmann" / "bachmann_dd_from_split_tables_compat.cdf"),
        help="Bachmann baseline CDF/CDL file.",
    )
    parser.add_argument(
        "--krstic-total-cdf",
        default=str(krstic_dir / "krstic_dd_total_elastic_integral_priority.cdf"),
        help="Krstic total-elastic CDF/CDL file.",
    )
    parser.add_argument(
        "--krstic-pure-cdf",
        default=str(krstic_dir / "krstic_dd_pure_dcs_compat.cdf"),
        help="Krstic pure-elastic CDF/CDL file.",
    )
    parser.add_argument(
        "--memo",
        default=str(this_dir / "DpD_fit_memo_v2.md"),
        help="Markdown memo holding the manually curated Krstic DCS coefficients.",
    )
    parser.add_argument(
        "--coeff-json",
        default=str(krstic_dir / "krstic_dd_dcs_coeffs.json"),
        help="Shared Krstic DCS coefficient JSON.",
    )
    parser.add_argument(
        "--metrics-csv",
        default=str(krstic_dir / "krstic_bachmann_angle_debug_metrics.csv"),
        help="Per-energy comparison metrics CSV output path.",
    )
    parser.add_argument(
        "--summary-json",
        default=str(krstic_dir / "krstic_bachmann_angle_debug_summary.json"),
        help="Summary JSON output path.",
    )
    parser.add_argument(
        "--samples-figure",
        default=str(this_dir / "figure" / "krstic_bachmann_angle_debug_samples.png"),
        help="Representative angle-curve comparison figure.",
    )
    parser.add_argument(
        "--native-metrics-csv",
        default=str(krstic_dir / "krstic_bachmann_native_dcs_angle_metrics.csv"),
        help="Native DCS-energy comparison metrics CSV output path.",
    )
    parser.add_argument(
        "--native-samples-figure",
        default=str(this_dir / "figure" / "krstic_bachmann_native_dcs_angle_samples.png"),
        help="Representative native DCS-energy angle-curve comparison figure.",
    )
    parser.add_argument(
        "--transport-figure",
        default=str(this_dir / "figure" / "krstic_bachmann_angle_debug_transport.png"),
        help="Transport-ratio comparison figure.",
    )
    parser.add_argument(
        "--sample-energies",
        type=float,
        nargs="*",
        default=[
            9.931816931229604e-04,
            1.9873367043051808,
            9.932710794753412,
            19.8644284078137,
            49.65958124699457,
            99.31816931229602,
        ],
        help="Requested sample energies on the runtime lab-energy axis.",
    )
    parser.add_argument(
        "--negative-pure-policy",
        choices=("assert", "warn-clip", "clip"),
        default="warn-clip",
        help="How to handle locally negative pure-elastic kernels during direct rebuild.",
    )
    parser.add_argument(
        "--theta-log-points",
        type=int,
        default=8192,
        help="Log-spaced theta points below the split angle for the pure-elastic rebuild.",
    )
    parser.add_argument(
        "--theta-linear-points",
        type=int,
        default=8192,
        help="Linearly spaced theta points above the split angle for the pure-elastic rebuild.",
    )
    parser.add_argument(
        "--theta-points-total",
        type=int,
        default=16384,
        help="Theta grid points for the total-elastic direct inverse CDF rebuild.",
    )
    return parser.parse_args()


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"No rows to write for {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def to_cdf_compatible_angle_order(angle_grid: np.ndarray) -> np.ndarray:
    return np.asarray(angle_grid, dtype=float)[:, ::-1]


def select_sample_indices(runtime_energy: np.ndarray, requested: list[float]) -> list[int]:
    selected: list[int] = []
    for target in requested:
        index = int(np.argmin(np.abs(runtime_energy - float(target))))
        if index not in selected:
            selected.append(index)
    return selected


def angle_at_probability(theta_row: np.ndarray, prob_axis: np.ndarray, probability: float) -> float:
    return float(np.interp(float(probability), prob_axis, theta_row))


def rms_angle_difference(lhs: np.ndarray, rhs: np.ndarray) -> float:
    diff = np.asarray(lhs, dtype=float) - np.asarray(rhs, dtype=float)
    return float(np.sqrt(np.mean(diff**2)))


def max_abs_angle_difference(lhs: np.ndarray, rhs: np.ndarray) -> float:
    return float(np.max(np.abs(np.asarray(lhs, dtype=float) - np.asarray(rhs, dtype=float))))


def worst_row(rows: list[dict[str, object]], key: str) -> dict[str, object]:
    return max(rows, key=lambda row: abs(float(row[key])))


def interpolate_angle_grid_by_energy(
    energy_in: np.ndarray,
    angle_grid_in: np.ndarray,
    energy_out: np.ndarray,
) -> np.ndarray:
    energy_src = np.asarray(energy_in, dtype=float)
    angle_src = np.asarray(angle_grid_in, dtype=float)
    energy_dst = np.asarray(energy_out, dtype=float)
    if angle_src.shape[0] != energy_src.size:
        raise ValueError("Angle-grid row count must match the energy axis length.")

    log_energy_src = np.log(energy_src)
    log_energy_dst = np.log(np.clip(energy_dst, energy_src[0], energy_src[-1]))
    angle_out = np.zeros((energy_dst.size, angle_src.shape[1]), dtype=float)
    for prob_index in range(angle_src.shape[1]):
        angle_out[:, prob_index] = np.interp(
            log_energy_dst,
            log_energy_src,
            angle_src[:, prob_index],
            left=float(angle_src[0, prob_index]),
            right=float(angle_src[-1, prob_index]),
        )
    return angle_out


def evaluate_transport_ratio_grid(prob_axis: np.ndarray, angle_grid: np.ndarray) -> np.ndarray:
    prob_values = np.asarray(prob_axis, dtype=float)
    return np.array(
        [
            float(np.trapz(1.0 - np.cos(np.asarray(theta_row, dtype=float)), x=prob_values))
            for theta_row in np.asarray(angle_grid, dtype=float)
        ],
        dtype=float,
    )


def plot_angle_samples(
    runtime_energy_lab: np.ndarray,
    prob_axis: np.ndarray,
    bachmann_angle: np.ndarray,
    krstic_total_angle: np.ndarray,
    krstic_pure_angle: np.ndarray,
    sample_indices: list[int],
    output_path: Path,
    *,
    bachmann_label: str,
    krstic_total_label: str,
    krstic_pure_label: str,
    figure_title: str,
) -> None:
    n_panels = len(sample_indices)
    if n_panels == 6:
        ncols = 3
        panel_width = 5.2
        panel_height = 4.2
        axes_box_aspect = 0.8
    else:
        ncols = min(3, max(1, n_panels))
        panel_width = 5.2
        panel_height = 3.6
        axes_box_aspect = None
    nrows = int(np.ceil(n_panels / ncols))
    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(panel_width * ncols, panel_height * nrows),
    )
    axes_array = np.atleast_1d(axes).ravel()
    probability = 1.0 - prob_axis

    for axis, energy_index in zip(axes_array, sample_indices, strict=False):
        energy_lab = float(runtime_energy_lab[energy_index])
        axis.plot(
            bachmann_angle[energy_index, ::-1],
            probability[::-1],
            lw=2.0,
            label=bachmann_label,
        )
        axis.plot(
            krstic_total_angle[energy_index, ::-1],
            probability[::-1],
            lw=1.8,
            label=krstic_total_label,
        )
        axis.plot(
            krstic_pure_angle[energy_index, ::-1],
            probability[::-1],
            lw=1.8,
            label=krstic_pure_label,
        )
        axis.set_title(f"E = {energy_lab:.4g} eV/amu")
        axis.set_xlim(0.0, np.pi)
        axis.set_ylim(0.0, 1.0)
        axis.set_xlabel("Scattering angle [rad]")
        axis.set_ylabel("Probability")
        if axes_box_aspect is not None:
            axis.set_box_aspect(axes_box_aspect)
        axis.grid(True, alpha=0.3)

    for axis in axes_array[n_panels:]:
        axis.axis("off")

    handles, labels = axes_array[0].get_legend_handles_labels()
    fig.suptitle(figure_title, y=0.99)
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.95),
        ncol=3,
        frameon=False,
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.86))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def plot_transport_ratio(
    runtime_energy_lab: np.ndarray,
    bachmann_ratio: np.ndarray,
    krstic_total_ratio: np.ndarray,
    krstic_pure_ratio: np.ndarray,
    pure_fit_range_lab: tuple[float, float],
    output_path: Path,
) -> None:
    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    ax.plot(runtime_energy_lab, bachmann_ratio, lw=2.0, label="Bachmann")
    ax.plot(runtime_energy_lab, krstic_total_ratio, lw=1.8, label="Krstic total")
    ax.plot(runtime_energy_lab, krstic_pure_ratio, lw=1.8, label="Krstic pure")
    ax.axvline(pure_fit_range_lab[0], color="0.4", lw=1.0, ls="--", label="pure DCS fit lower bound")
    ax.set_xscale("log")
    ax.set_xlabel("E [eV/amu]")
    ax.set_ylabel(r"$\int (1-\cos\theta)\,dR$")
    ax.set_title("Transport ratio from scattering-angle tables")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=9)
    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    this_dir = Path(__file__).resolve().parent
    krstic_dir = this_dir / "Krstic"
    sys.path.insert(0, str(this_dir))

    from cdf_compat import (
        ElasticCdf,
        angle_axes,
        evaluate_angle_transport_ratio,
        reshape_scattering_angles,
    )
    from krstic_dcs import (
        build_pure_dcs_first_observables_by_energy,
        build_tabulated_observables,
        load_or_build_coefficients,
    )

    krstic_total_module = load_path_module(
        "krstic_total_debug_module",
        krstic_dir / "krstic_dcs.py",
    )

    bachmann_cdf = ElasticCdf.from_file(args.bachmann_cdf)
    krstic_total_cdf = ElasticCdf.from_file(args.krstic_total_cdf)
    krstic_pure_cdf = ElasticCdf.from_file(args.krstic_pure_cdf)

    prob_axis, runtime_energy_lab = angle_axes(bachmann_cdf)
    total_prob_axis, total_runtime_energy_lab = angle_axes(krstic_total_cdf)
    pure_prob_axis, pure_runtime_energy_lab = angle_axes(krstic_pure_cdf)
    if not np.allclose(prob_axis, total_prob_axis) or not np.allclose(prob_axis, pure_prob_axis):
        raise ValueError("Probability axes differ between Bachmann and Krstic tables.")
    if not np.allclose(runtime_energy_lab, total_runtime_energy_lab) or not np.allclose(
        runtime_energy_lab,
        pure_runtime_energy_lab,
    ):
        raise ValueError("Runtime energy axes differ between Bachmann and Krstic tables.")

    bachmann_angle = reshape_scattering_angles(bachmann_cdf.get_block(5), bachmann_cdf)
    krstic_total_angle = reshape_scattering_angles(krstic_total_cdf.get_block(5), krstic_total_cdf)
    krstic_pure_angle = reshape_scattering_angles(krstic_pure_cdf.get_block(5), krstic_pure_cdf)

    bachmann_ratio = evaluate_angle_transport_ratio(bachmann_cdf.get_block(5), bachmann_cdf)
    krstic_total_ratio = evaluate_angle_transport_ratio(
        krstic_total_cdf.get_block(5),
        krstic_total_cdf,
    )
    krstic_pure_ratio = evaluate_angle_transport_ratio(
        krstic_pure_cdf.get_block(5),
        krstic_pure_cdf,
    )

    runtime_energy_cm = 0.5 * runtime_energy_lab

    total_fits = krstic_total_module.parse_krstic_dcs_markdown(args.memo)["Elastic"]
    total_fit_energy_cm, total_inverse_native, _ = krstic_total_module.build_tabulated_inverse_cdf(
        total_fits,
        prob_axis,
        theta_points=args.theta_points_total,
    )
    total_runtime_inverse, _ = krstic_total_module.build_dcs_first_inverse_cdf_by_energy(
        total_fits,
        prob_axis,
        runtime_energy_cm,
        theta_points=args.theta_points_total,
    )
    total_runtime_rebuilt = to_cdf_compatible_angle_order(
        total_runtime_inverse
    )

    energy_slices = load_or_build_coefficients(
        coeff_json_path=args.coeff_json,
        markdown_path=args.memo,
    )
    pure_tabulated = build_tabulated_observables(
        energy_slices,
        prob_axis,
        log_points=args.theta_log_points,
        linear_points=args.theta_linear_points,
        negative_policy=args.negative_pure_policy,
    )
    pure_runtime_dcs = build_pure_dcs_first_observables_by_energy(
        energy_slices,
        prob_axis,
        runtime_energy_cm,
        log_points=args.theta_log_points,
        linear_points=args.theta_linear_points,
        negative_policy=args.negative_pure_policy,
    )
    pure_runtime_rebuilt = to_cdf_compatible_angle_order(
        np.asarray(pure_runtime_dcs["inverse_cdf_rad"], dtype=float)
    )
    if not np.allclose(total_fit_energy_cm, pure_tabulated["energy_cm_ev"]):
        raise ValueError("Total and pure Krstic native DCS energy grids differ.")
    native_energy_lab = 2.0 * np.asarray(total_fit_energy_cm, dtype=float)
    total_native_angle = to_cdf_compatible_angle_order(total_inverse_native)
    pure_native_angle = to_cdf_compatible_angle_order(
        np.asarray(pure_tabulated["inverse_cdf_rad"], dtype=float)
    )
    total_native_ratio = evaluate_transport_ratio_grid(prob_axis, total_native_angle)
    pure_native_ratio = evaluate_transport_ratio_grid(prob_axis, pure_native_angle)
    native_overlap_mask = native_energy_lab <= float(runtime_energy_lab[-1])
    native_overlap_energy_lab = native_energy_lab[native_overlap_mask]
    total_native_angle_overlap = total_native_angle[native_overlap_mask]
    pure_native_angle_overlap = pure_native_angle[native_overlap_mask]
    total_native_ratio_overlap = total_native_ratio[native_overlap_mask]
    pure_native_ratio_overlap = pure_native_ratio[native_overlap_mask]
    bachmann_native_angle_overlap = interpolate_angle_grid_by_energy(
        runtime_energy_lab,
        bachmann_angle,
        native_overlap_energy_lab,
    )
    bachmann_native_ratio_overlap = evaluate_transport_ratio_grid(
        prob_axis,
        bachmann_native_angle_overlap,
    )
    native_plot_sample_indices = select_sample_indices(native_energy_lab, args.sample_energies)
    bachmann_native_plot_angle = interpolate_angle_grid_by_energy(
        runtime_energy_lab,
        bachmann_angle,
        native_energy_lab,
    )

    pure_fit_range_lab = (
        2.0 * float(pure_tabulated["energy_cm_ev"][0]),
        2.0 * float(pure_tabulated["energy_cm_ev"][-1]),
    )
    sample_indices = select_sample_indices(runtime_energy_lab, args.sample_energies)

    rows: list[dict[str, object]] = []
    for index, energy_lab in enumerate(runtime_energy_lab):
        row = {
            "energy_lab_ev": float(energy_lab),
            "energy_cm_ev": float(0.5 * energy_lab),
            "pure_fit_energy_below_lower_bound": bool(energy_lab < pure_fit_range_lab[0]),
            "bachmann_transport_ratio": float(bachmann_ratio[index]),
            "krstic_total_transport_ratio": float(krstic_total_ratio[index]),
            "krstic_pure_transport_ratio": float(krstic_pure_ratio[index]),
            "total_minus_bachmann_transport_ratio": float(
                krstic_total_ratio[index] - bachmann_ratio[index]
            ),
            "pure_minus_bachmann_transport_ratio": float(
                krstic_pure_ratio[index] - bachmann_ratio[index]
            ),
            "total_vs_bachmann_rms_angle_diff_rad": rms_angle_difference(
                krstic_total_angle[index],
                bachmann_angle[index],
            ),
            "pure_vs_bachmann_rms_angle_diff_rad": rms_angle_difference(
                krstic_pure_angle[index],
                bachmann_angle[index],
            ),
            "total_vs_bachmann_max_abs_angle_diff_rad": max_abs_angle_difference(
                krstic_total_angle[index],
                bachmann_angle[index],
            ),
            "pure_vs_bachmann_max_abs_angle_diff_rad": max_abs_angle_difference(
                krstic_pure_angle[index],
                bachmann_angle[index],
            ),
            "total_stored_vs_rebuilt_max_abs_angle_diff_rad": max_abs_angle_difference(
                krstic_total_angle[index],
                total_runtime_rebuilt[index],
            ),
            "pure_stored_vs_rebuilt_max_abs_angle_diff_rad": max_abs_angle_difference(
                krstic_pure_angle[index],
                pure_runtime_rebuilt[index],
            ),
        }
        for probability in (0.1, 0.2, 0.5, 0.9):
            label = str(probability).replace(".", "p")
            row[f"bachmann_theta_at_{label}_rad"] = angle_at_probability(
                bachmann_angle[index],
                prob_axis,
                probability,
            )
            row[f"krstic_total_theta_at_{label}_rad"] = angle_at_probability(
                krstic_total_angle[index],
                prob_axis,
                probability,
            )
            row[f"krstic_pure_theta_at_{label}_rad"] = angle_at_probability(
                krstic_pure_angle[index],
                prob_axis,
                probability,
            )
        rows.append(row)

    native_rows: list[dict[str, object]] = []
    for index, energy_lab in enumerate(native_overlap_energy_lab):
        row = {
            "energy_lab_ev": float(energy_lab),
            "energy_cm_ev": float(0.5 * energy_lab),
            "bachmann_transport_ratio": float(bachmann_native_ratio_overlap[index]),
            "krstic_total_native_transport_ratio": float(total_native_ratio_overlap[index]),
            "krstic_pure_native_transport_ratio": float(pure_native_ratio_overlap[index]),
            "total_native_minus_bachmann_transport_ratio": float(
                total_native_ratio_overlap[index] - bachmann_native_ratio_overlap[index]
            ),
            "pure_native_minus_bachmann_transport_ratio": float(
                pure_native_ratio_overlap[index] - bachmann_native_ratio_overlap[index]
            ),
            "total_native_vs_bachmann_rms_angle_diff_rad": rms_angle_difference(
                total_native_angle_overlap[index],
                bachmann_native_angle_overlap[index],
            ),
            "pure_native_vs_bachmann_rms_angle_diff_rad": rms_angle_difference(
                pure_native_angle_overlap[index],
                bachmann_native_angle_overlap[index],
            ),
            "total_native_vs_bachmann_max_abs_angle_diff_rad": max_abs_angle_difference(
                total_native_angle_overlap[index],
                bachmann_native_angle_overlap[index],
            ),
            "pure_native_vs_bachmann_max_abs_angle_diff_rad": max_abs_angle_difference(
                pure_native_angle_overlap[index],
                bachmann_native_angle_overlap[index],
            ),
        }
        for probability in (0.1, 0.2, 0.5, 0.9):
            label = str(probability).replace(".", "p")
            row[f"bachmann_theta_at_{label}_rad"] = angle_at_probability(
                bachmann_native_angle_overlap[index],
                prob_axis,
                probability,
            )
            row[f"krstic_total_native_theta_at_{label}_rad"] = angle_at_probability(
                total_native_angle_overlap[index],
                prob_axis,
                probability,
            )
            row[f"krstic_pure_native_theta_at_{label}_rad"] = angle_at_probability(
                pure_native_angle_overlap[index],
                prob_axis,
                probability,
            )
        native_rows.append(row)

    metrics_path = Path(args.metrics_csv)
    native_metrics_path = Path(args.native_metrics_csv)
    summary_path = Path(args.summary_json)
    samples_figure_path = Path(args.samples_figure)
    native_samples_figure_path = Path(args.native_samples_figure)
    transport_figure_path = Path(args.transport_figure)
    for output_path in (
        metrics_path,
        native_metrics_path,
        summary_path,
        samples_figure_path,
        native_samples_figure_path,
        transport_figure_path,
    ):
        output_path.parent.mkdir(parents=True, exist_ok=True)

    write_csv(metrics_path, rows)
    write_csv(native_metrics_path, native_rows)

    summary = {
        "bachmann_cdf": str(Path(args.bachmann_cdf).resolve()),
        "krstic_total_cdf": str(Path(args.krstic_total_cdf).resolve()),
        "krstic_pure_cdf": str(Path(args.krstic_pure_cdf).resolve()),
        "memo": str(Path(args.memo).resolve()),
        "coeff_json": str(Path(args.coeff_json).resolve()),
        "runtime_energy_lab_ev_range": [
            float(runtime_energy_lab[0]),
            float(runtime_energy_lab[-1]),
        ],
        "native_overlap_energy_lab_ev_range": [
            float(native_overlap_energy_lab[0]),
            float(native_overlap_energy_lab[-1]),
        ],
        "native_overlap_row_count": int(native_overlap_energy_lab.size),
        "pure_fit_energy_lab_ev_range": list(pure_fit_range_lab),
        "runtime_energy_below_pure_fit_range_count": int(
            np.count_nonzero(runtime_energy_lab < pure_fit_range_lab[0])
        ),
        "sample_energies_lab_ev": [float(runtime_energy_lab[index]) for index in sample_indices],
        "native_sample_energies_lab_ev": [
            float(native_energy_lab[index]) for index in native_plot_sample_indices
        ],
        "max_total_vs_bachmann_rms_angle_diff_rad": float(
            max(float(row["total_vs_bachmann_rms_angle_diff_rad"]) for row in rows)
        ),
        "max_pure_vs_bachmann_rms_angle_diff_rad": float(
            max(float(row["pure_vs_bachmann_rms_angle_diff_rad"]) for row in rows)
        ),
        "max_total_native_vs_bachmann_rms_angle_diff_rad": float(
            max(float(row["total_native_vs_bachmann_rms_angle_diff_rad"]) for row in native_rows)
        ),
        "max_pure_native_vs_bachmann_rms_angle_diff_rad": float(
            max(float(row["pure_native_vs_bachmann_rms_angle_diff_rad"]) for row in native_rows)
        ),
        "max_total_vs_bachmann_transport_ratio_delta": float(
            max(abs(float(row["total_minus_bachmann_transport_ratio"])) for row in rows)
        ),
        "max_pure_vs_bachmann_transport_ratio_delta": float(
            max(abs(float(row["pure_minus_bachmann_transport_ratio"])) for row in rows)
        ),
        "max_total_native_vs_bachmann_transport_ratio_delta": float(
            max(
                abs(float(row["total_native_minus_bachmann_transport_ratio"]))
                for row in native_rows
            )
        ),
        "max_pure_native_vs_bachmann_transport_ratio_delta": float(
            max(
                abs(float(row["pure_native_minus_bachmann_transport_ratio"]))
                for row in native_rows
            )
        ),
        "max_total_stored_vs_rebuilt_angle_diff_rad": float(
            max(float(row["total_stored_vs_rebuilt_max_abs_angle_diff_rad"]) for row in rows)
        ),
        "max_pure_stored_vs_rebuilt_angle_diff_rad": float(
            max(float(row["pure_stored_vs_rebuilt_max_abs_angle_diff_rad"]) for row in rows)
        ),
        "worst_total_vs_bachmann_rms_row": worst_row(rows, "total_vs_bachmann_rms_angle_diff_rad"),
        "worst_pure_vs_bachmann_rms_row": worst_row(rows, "pure_vs_bachmann_rms_angle_diff_rad"),
        "worst_total_native_vs_bachmann_rms_row": worst_row(
            native_rows,
            "total_native_vs_bachmann_rms_angle_diff_rad",
        ),
        "worst_pure_native_vs_bachmann_rms_row": worst_row(
            native_rows,
            "pure_native_vs_bachmann_rms_angle_diff_rad",
        ),
        "worst_total_transport_ratio_delta_row": worst_row(
            rows,
            "total_minus_bachmann_transport_ratio",
        ),
        "worst_pure_transport_ratio_delta_row": worst_row(
            rows,
            "pure_minus_bachmann_transport_ratio",
        ),
        "worst_total_native_transport_ratio_delta_row": worst_row(
            native_rows,
            "total_native_minus_bachmann_transport_ratio",
        ),
        "worst_pure_native_transport_ratio_delta_row": worst_row(
            native_rows,
            "pure_native_minus_bachmann_transport_ratio",
        ),
    }
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    plot_angle_samples(
        runtime_energy_lab,
        prob_axis,
        bachmann_angle,
        krstic_total_angle,
        krstic_pure_angle,
        sample_indices,
        samples_figure_path,
        bachmann_label="Bachmann table",
        krstic_total_label="Krstic total-elastic table",
        krstic_pure_label="Krstic pure-elastic table",
        figure_title="Bachmann vs Krstic runtime scattering-angle CDF samples",
    )
    plot_angle_samples(
        native_energy_lab,
        prob_axis,
        bachmann_native_plot_angle,
        total_native_angle,
        pure_native_angle,
        native_plot_sample_indices,
        native_samples_figure_path,
        bachmann_label="Bachmann table (interpolated/clamped)",
        krstic_total_label="Krstic total-elastic native DCS",
        krstic_pure_label="Krstic pure-elastic native DCS",
        figure_title="Bachmann vs Krstic native-DCS scattering-angle CDF samples",
    )
    plot_transport_ratio(
        runtime_energy_lab,
        bachmann_ratio,
        krstic_total_ratio,
        krstic_pure_ratio,
        pure_fit_range_lab,
        transport_figure_path,
    )

    print(f"Wrote {metrics_path}")
    print(f"Wrote {native_metrics_path}")
    print(f"Wrote {summary_path}")
    print(f"Wrote {samples_figure_path}")
    print(f"Wrote {native_samples_figure_path}")
    print(f"Wrote {transport_figure_path}")
    print(f"  pure fit range [lab eV]        : {pure_fit_range_lab[0]:.6g} .. {pure_fit_range_lab[1]:.6g}")
    print(
        "  runtime rows below pure fit    : "
        f"{int(np.count_nonzero(runtime_energy_lab < pure_fit_range_lab[0]))}"
    )
    print(
        "  native overlap [lab eV]        : "
        f"{native_overlap_energy_lab[0]:.6g} .. {native_overlap_energy_lab[-1]:.6g}"
    )
    print(
        "  max stored-vs-rebuilt diff [rad]: "
        f"total={summary['max_total_stored_vs_rebuilt_angle_diff_rad']:.3e}, "
        f"pure={summary['max_pure_stored_vs_rebuilt_angle_diff_rad']:.3e}"
    )
    print(
        "  max native-vs-Bachmann RMS [rad]: "
        f"total={summary['max_total_native_vs_bachmann_rms_angle_diff_rad']:.3e}, "
        f"pure={summary['max_pure_native_vs_bachmann_rms_angle_diff_rad']:.3e}"
    )


if __name__ == "__main__":
    main()
