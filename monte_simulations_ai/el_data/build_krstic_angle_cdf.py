#!/usr/bin/env python3
"""Krstic full-rebuild workflow: rebuild a dd_00_elastic-compatible CDF from direct DCS fits."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"No rows to write for {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def to_cdf_compatible_angle_order(angle_grid: np.ndarray) -> np.ndarray:
    """Match the legacy dd_00_elastic.cdf convention: pi -> 0 along probability."""
    return np.asarray(angle_grid, dtype=float)[:, ::-1]


def parse_args() -> argparse.Namespace:
    this_dir = Path(__file__).resolve().parent
    out_dir = this_dir / "Krstic"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--template-cdf",
        default=str(this_dir / "dd_00_elastic.cdf"),
        help="Template dd_00_elastic.cdf/CDL file.",
    )
    parser.add_argument(
        "--memo",
        default=str(this_dir / "DpD_fit_memo_v2.md"),
        help="Markdown memo that holds the manually curated Krstic DCS coefficients.",
    )
    parser.add_argument(
        "--coeff-json",
        default=str(out_dir / "krstic_dd_dcs_coeffs.json"),
        help="Shared output/input Krstic DCS coefficient JSON.",
    )
    parser.add_argument(
        "--output-cdf",
        default=str(out_dir / "krstic_dd_pure_dcs_compat.cdf"),
        help="CDF/CDL output path.",
    )
    parser.add_argument(
        "--validation-json",
        default=str(out_dir / "krstic_pure_dcs_validation.json"),
        help="Validation summary JSON output path.",
    )
    parser.add_argument(
        "--validation-csv",
        default=str(out_dir / "krstic_pure_dcs_validation.csv"),
        help="Per-energy validation CSV output path.",
    )
    parser.add_argument(
        "--angle-csv",
        default=str(out_dir / "krstic_pure_scattering_angle_compat.csv"),
        help="CSV for the 251 x 51 angle table.",
    )
    parser.add_argument(
        "--transport-csv",
        default=str(out_dir / "krstic_pure_transport_from_dcs.csv"),
        help="CSV for reaction_rate and I_1_x tables.",
    )
    parser.add_argument(
        "--negative-pure-policy",
        choices=("assert", "warn-clip", "clip"),
        default="warn-clip",
        help="How to handle locally negative pure-elastic kernels.",
    )
    parser.add_argument(
        "--xi-points",
        type=int,
        default=1000,
        help="Quadrature points for the xi integration in the I-kernel tables.",
    )
    parser.add_argument(
        "--theta-log-points",
        type=int,
        default=8192,
        help="Log-spaced theta points below the split angle.",
    )
    parser.add_argument(
        "--theta-linear-points",
        type=int,
        default=8192,
        help="Linearly spaced theta points above the split angle.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    this_dir = Path(__file__).resolve().parent
    sys.path.insert(0, str(this_dir))

    from cdf_compat import (
        ElasticCdf,
        angle_axes,
        apply_blocks,
        compute_transport_tables_from_sigma_momentum,
        cross_section_axis,
        evaluate_angle_transport_ratio,
        transport_axes,
    )
    from krstic_dcs import (
        BOHR_AREA_CM2,
        build_pure_dcs_first_observables_by_energy,
        build_tabulated_observables,
        load_or_build_coefficients,
        write_coefficients_json,
    )
    sys.path.append(str(this_dir / "Krstic"))
    from generate_krstic_integral_data import MODELS

    template = ElasticCdf.from_file(args.template_cdf)
    coeff_json_path = Path(args.coeff_json)
    memo_path = Path(args.memo)
    for output_path in (
        coeff_json_path,
        Path(args.output_cdf),
        Path(args.validation_json),
        Path(args.validation_csv),
        Path(args.angle_csv),
        Path(args.transport_csv),
    ):
        output_path.parent.mkdir(parents=True, exist_ok=True)
    write_coefficients_json(memo_path, coeff_json_path)
    energy_slices = load_or_build_coefficients(
        coeff_json_path=coeff_json_path,
        markdown_path=memo_path,
    )

    prob_axis, angle_energy_lab = angle_axes(template)
    tabulated = build_tabulated_observables(
        energy_slices,
        prob_axis,
        log_points=args.theta_log_points,
        linear_points=args.theta_linear_points,
        negative_policy=args.negative_pure_policy,
    )

    sigma_energy_lab = cross_section_axis(template)
    sigma_energy_cm = 0.5 * sigma_energy_lab
    angle_energy_cm = 0.5 * angle_energy_lab

    sigma_runtime = build_pure_dcs_first_observables_by_energy(
        energy_slices,
        prob_axis,
        sigma_energy_cm,
        log_points=args.theta_log_points,
        linear_points=args.theta_linear_points,
        negative_policy=args.negative_pure_policy,
    )
    angle_runtime = build_pure_dcs_first_observables_by_energy(
        energy_slices,
        prob_axis,
        angle_energy_cm,
        log_points=args.theta_log_points,
        linear_points=args.theta_linear_points,
        negative_policy=args.negative_pure_policy,
    )
    sigma_t_pure_cm2 = np.asarray(sigma_runtime["sigma_t_pure_cm2"], dtype=float)
    sigma_mt_pure_cm2 = np.asarray(sigma_runtime["sigma_mt_pure_cm2"], dtype=float)
    angle_table = np.asarray(angle_runtime["inverse_cdf_rad"], dtype=float)
    angle_table_compat = to_cdf_compatible_angle_order(angle_table)
    scattering_angle = angle_table_compat.ravel()

    transport = compute_transport_tables_from_sigma_momentum(
        template,
        cross_section=sigma_t_pure_cm2,
        sigma_momentum=sigma_mt_pure_cm2,
        sigma_energy=sigma_energy_lab,
        xi_points=args.xi_points,
    )

    apply_blocks(
        template,
        cross_section=sigma_t_pure_cm2,
        scattering_angle=scattering_angle,
        reaction_rate=transport["reaction_rate"],
        i_1_0=transport["I_1_0"],
        i_1_1_up=transport["I_1_1_up"],
        i_1_2_up2=transport["I_1_2_up2"],
        sigv_max=float(transport["sigv_max"]),
        angle_min=0.0,
    )
    data_version = (
        "Data Version 3.4 D + D^+ pure-elastic cross section and angle CDF "
        "rebuilt DCS-first from Krstic direct DCS fits."
    )
    template.write(
        args.output_cdf,
        xs_name="dd_00_el_krstic_pure_dcs",
        description=(
            "Krstic D+ + D pure-elastic tables rebuilt DCS-first from direct DCS fits."
        ),
        data_version=data_version,
    )

    validation_rows: list[dict[str, object]] = []
    for row in tabulated["rows"]:
        energy_cm = float(row["energy_cm_ev"])
        integral_sigma_t = MODELS["elastic"].evaluate_au(energy_cm) * BOHR_AREA_CM2
        integral_sigma_mt = MODELS["momentum_transfer"].evaluate_au(energy_cm) * BOHR_AREA_CM2
        integral_sigma_vi = MODELS["viscosity"].evaluate_au(energy_cm) * BOHR_AREA_CM2
        integral_sigma_se = MODELS["spin_exchange"].evaluate_au(energy_cm) * BOHR_AREA_CM2
        pure_reference_sigma_t = integral_sigma_t - integral_sigma_se
        pure_ratio = float(row["sigma_mt_pure_cm2"]) / float(row["sigma_t_pure_cm2"])
        validation_rows.append(
            {
                "energy_cm_ev": energy_cm,
                "energy_lab_ev": 2.0 * energy_cm,
                "sigma_t_elastic_total_dcs_cm2": float(row["sigma_t_elastic_total_cm2"]),
                "sigma_t_elastic_integral_cm2": float(integral_sigma_t),
                "sigma_t_elastic_relative_error": float(
                    row["sigma_t_elastic_total_cm2"] / integral_sigma_t - 1.0
                ),
                "sigma_mt_elastic_total_dcs_cm2": float(row["sigma_mt_elastic_total_cm2"]),
                "sigma_mt_integral_cm2": float(integral_sigma_mt),
                "sigma_mt_elastic_relative_error": float(
                    row["sigma_mt_elastic_total_cm2"] / integral_sigma_mt - 1.0
                ),
                "sigma_vi_elastic_total_dcs_cm2": float(row["sigma_vi_elastic_total_cm2"]),
                "sigma_vi_integral_cm2": float(integral_sigma_vi),
                "sigma_vi_elastic_relative_error": float(
                    row["sigma_vi_elastic_total_cm2"] / integral_sigma_vi - 1.0
                ),
                "sigma_t_pure_dcs_cm2": float(row["sigma_t_pure_cm2"]),
                "sigma_t_elastic_minus_se_cm2": float(pure_reference_sigma_t),
                "sigma_t_pure_relative_error_vs_el_minus_se": float(
                    row["sigma_t_pure_cm2"] / pure_reference_sigma_t - 1.0
                ),
                "sigma_mt_pure_dcs_cm2": float(row["sigma_mt_pure_cm2"]),
                "sigma_vi_pure_dcs_cm2": float(row["sigma_vi_pure_cm2"]),
                "cdf_transport_ratio": float(row["transport_ratio_from_cdf"]),
                "moment_transport_ratio": pure_ratio,
                "cdf_vs_moment_ratio_error": float(
                    row["transport_ratio_from_cdf"] / pure_ratio - 1.0
                ),
                "theta_min_used_rad": float(row["theta_min_used_rad"]),
                "pole_count": len(row["pole_thetas_rad"]),
                "negative_min_raw_au": float(row["negative_min_raw_au"]),
                "clipped_negative_points": int(row["clipped_negative_points"]),
            }
        )
    write_csv(Path(args.validation_csv), validation_rows)

    runtime_angle_rows: list[dict[str, object]] = []
    for energy_index, energy_lab in enumerate(angle_energy_lab):
        for prob_index, probability in enumerate(prob_axis):
            runtime_angle_rows.append(
                {
                    "energy_lab_ev": float(energy_lab),
                    "energy_cm_ev": float(0.5 * energy_lab),
                    "probability": float(probability),
                    "theta_rad": float(angle_table_compat[energy_index, prob_index]),
                }
            )
    write_csv(Path(args.angle_csv), runtime_angle_rows)

    transport_rows: list[dict[str, object]] = []
    energy_transport_lab, temp_transport = transport_axes(template)
    rr_table = transport["reaction_rate"].reshape(temp_transport.size, energy_transport_lab.size)
    i10_table = transport["I_1_0"].reshape(temp_transport.size, energy_transport_lab.size)
    i11_table = transport["I_1_1_up"].reshape(temp_transport.size, energy_transport_lab.size)
    i12_table = transport["I_1_2_up2"].reshape(temp_transport.size, energy_transport_lab.size)
    for temp_index, temp_lab in enumerate(temp_transport):
        for energy_index, energy_lab in enumerate(energy_transport_lab):
            transport_rows.append(
                {
                    "projectile_energy_lab_ev": float(energy_lab),
                    "projectile_energy_cm_ev": float(0.5 * energy_lab),
                    "background_temperature_lab_ev": float(temp_lab),
                    "background_temperature_cm_ev": float(0.5 * temp_lab),
                    "reaction_rate_cm3_s": float(rr_table[temp_index, energy_index]),
                    "I_1_0_cm3_s": float(i10_table[temp_index, energy_index]),
                    "I_1_1_up_cm4_s2": float(i11_table[temp_index, energy_index]),
                    "I_1_2_up2_cm5_s3": float(i12_table[temp_index, energy_index]),
                }
            )
    write_csv(Path(args.transport_csv), transport_rows)

    runtime_angle_transport_ratio = evaluate_angle_transport_ratio(scattering_angle, template)
    runtime_angle_dcs_ratio = np.asarray(angle_runtime["transport_ratio_pure"], dtype=float)
    runtime_angle_cdf_ratio = np.asarray(angle_runtime["transport_ratio_from_cdf"], dtype=float)

    validation = {
        "coeff_json": str(coeff_json_path),
        "memo": str(memo_path),
        "data_version": data_version,
        "cross_section_source": "Krstic pure-elastic direct DCS",
        "sigma_momentum_source": "Krstic pure-elastic direct DCS",
        "angle_source": "Krstic pure-elastic direct DCS",
        "runtime_interpolation_policy": (
            "DCS-first: interpolate the elastic-total and spin-exchange "
            "p(theta,E)=2*pi*sin(theta)*d sigma/d Omega kernels separately in "
            "log(E)-log(p), subtract them to form pure elastic, clip negative "
            "pure points according to negative_pure_policy, then integrate to "
            "sigma_t/sigma_mt and inverse CDF."
        ),
        "below_tabulated_energy_policy": "Clamp to the lowest direct-DCS slice.",
        "negative_pure_policy": args.negative_pure_policy,
        "theta_log_points": args.theta_log_points,
        "theta_linear_points": args.theta_linear_points,
        "xi_points": args.xi_points,
        "fit_energy_range_cm_ev": [
            float(tabulated["energy_cm_ev"][0]),
            float(tabulated["energy_cm_ev"][-1]),
        ],
        "runtime_cross_section_range_lab_ev": [
            float(sigma_energy_lab[0]),
            float(sigma_energy_lab[-1]),
        ],
        "runtime_angle_range_lab_ev": [
            float(angle_energy_lab[0]),
            float(angle_energy_lab[-1]),
        ],
        "max_abs_sigma_t_elastic_relative_error": float(
            max(abs(row["sigma_t_elastic_relative_error"]) for row in validation_rows)
        ),
        "max_abs_sigma_mt_elastic_relative_error": float(
            max(abs(row["sigma_mt_elastic_relative_error"]) for row in validation_rows)
        ),
        "max_abs_sigma_vi_elastic_relative_error": float(
            max(abs(row["sigma_vi_elastic_relative_error"]) for row in validation_rows)
        ),
        "max_abs_sigma_t_pure_relative_error_vs_el_minus_se": float(
            max(abs(row["sigma_t_pure_relative_error_vs_el_minus_se"]) for row in validation_rows)
        ),
        "max_abs_cdf_vs_moment_ratio_error": float(
            max(abs(row["cdf_vs_moment_ratio_error"]) for row in validation_rows)
        ),
        "max_abs_runtime_angle_vs_dcs_moment_ratio_error": float(
            np.max(np.abs(runtime_angle_transport_ratio - runtime_angle_dcs_ratio))
        ),
        "max_abs_runtime_angle_vs_dcs_inverse_cdf_ratio_error": float(
            np.max(np.abs(runtime_angle_transport_ratio - runtime_angle_cdf_ratio))
        ),
        "runtime_cross_section_rows": sigma_runtime["rows"],
        "runtime_angle_rows": angle_runtime["rows"],
        "sigv_max_cm3_s": float(transport["sigv_max"]),
        "validation_rows": validation_rows,
    }
    Path(args.validation_json).write_text(
        json.dumps(validation, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    print(f"Wrote {coeff_json_path}")
    print(f"Wrote {args.validation_csv}")
    print(f"Wrote {args.validation_json}")
    print(f"Wrote {args.angle_csv}")
    print(f"Wrote {args.transport_csv}")
    print(f"Wrote {args.output_cdf}")
    print(f"  sigv_max: {float(transport['sigv_max']):.6e} cm^3/s")


if __name__ == "__main__":
    main()
