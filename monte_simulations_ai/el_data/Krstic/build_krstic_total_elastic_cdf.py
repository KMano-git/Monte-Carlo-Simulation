#!/usr/bin/env python3
"""Build a Krstic total-elastic dataset with integral-fit sigma and DCS-first angles."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np


THETA_FLOOR_RAD = 1.0e-8
THETA_POINTS = 16384


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


def evaluate_model_cm2(model, energy_cm_grid: np.ndarray, bohr_area_cm2: float) -> np.ndarray:
    return np.array(
        [model.evaluate_au(float(energy_cm)) * bohr_area_cm2 for energy_cm in energy_cm_grid],
        dtype=float,
    )


def serialize_shared_coefficients(
    fits_by_channel: dict[str, list[object]],
) -> dict[str, dict[str, dict[str, object]]]:
    elastic_by_energy = {float(fit.energy_ev): fit for fit in fits_by_channel["Elastic"]}
    spin_exchange_by_energy = {
        float(fit.energy_ev): fit for fit in fits_by_channel["Spin Exchange"]
    }
    energies = sorted(elastic_by_energy)
    if energies != sorted(spin_exchange_by_energy):
        raise ValueError("Elastic and spin-exchange DCS energies do not match.")

    serialized: dict[str, dict[str, dict[str, object]]] = {}
    for energy in energies:
        elastic_fit = elastic_by_energy[energy]
        spin_exchange_fit = spin_exchange_by_energy[energy]
        serialized[f"{energy:.4f}"] = {
            "elastic": {
                "a": list(elastic_fit.a_values),
                "b": list(elastic_fit.b_values),
                "A": float(elastic_fit.A),
                "B": float(elastic_fit.B),
                "C": float(elastic_fit.C),
            },
            "spin_exchange": {
                "a": list(spin_exchange_fit.a_values),
                "b": list(spin_exchange_fit.b_values),
                "A": float(spin_exchange_fit.A),
                "B": float(spin_exchange_fit.B),
                "C": float(spin_exchange_fit.C),
            },
        }
    return serialized


def parse_args() -> argparse.Namespace:
    this_dir = Path(__file__).resolve().parent
    el_data_dir = this_dir.parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--template-cdf",
        default=str(el_data_dir / "dd_00_elastic.cdf"),
        help="Template dd_00_elastic.cdf/CDL file.",
    )
    parser.add_argument(
        "--memo",
        default=str(el_data_dir / "DpD_fit_memo_v2.md"),
        help="Markdown memo that holds the manually curated Krstic DCS coefficients.",
    )
    parser.add_argument(
        "--output-cdf",
        default=str(this_dir / "krstic_dd_total_elastic_integral_priority.cdf"),
        help="Output CDF/CDL path.",
    )
    parser.add_argument(
        "--coeff-json",
        default=str(this_dir / "krstic_dd_dcs_coeffs.json"),
        help="Shared Krstic DCS coefficient JSON output path.",
    )
    parser.add_argument(
        "--validation-json",
        default=str(this_dir / "krstic_total_elastic_validation.json"),
        help="Validation summary JSON output path.",
    )
    parser.add_argument(
        "--validation-csv",
        default=str(this_dir / "krstic_total_elastic_validation.csv"),
        help="Per-energy validation CSV output path.",
    )
    parser.add_argument(
        "--angle-csv",
        default=str(this_dir / "krstic_total_elastic_scattering_angle_compat.csv"),
        help="CSV for the 251 x 51 angle table.",
    )
    parser.add_argument(
        "--transport-csv",
        default=str(this_dir / "krstic_total_elastic_transport.csv"),
        help="CSV for reaction_rate and I_1_x tables.",
    )
    parser.add_argument(
        "--xi-points",
        type=int,
        default=1000,
        help="Quadrature points for the xi integration in the transport tables.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    this_dir = Path(__file__).resolve().parent
    el_data_dir = this_dir.parent
    sys.path.insert(0, str(el_data_dir))
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
    from generate_krstic_integral_data import BOHR_AREA_CM2, MODELS
    from krstic_dcs import (
        build_dcs_first_inverse_cdf_by_energy,
        build_tabulated_inverse_cdf,
        parse_krstic_dcs_markdown,
    )

    template_path = Path(args.template_cdf)
    markdown_path = Path(args.memo)
    output_path = Path(args.output_cdf)
    coeff_json_out = Path(args.coeff_json)
    validation_json_out = Path(args.validation_json)
    validation_csv_out = Path(args.validation_csv)
    angle_csv_out = Path(args.angle_csv)
    transport_csv_out = Path(args.transport_csv)
    for output in (
        output_path,
        coeff_json_out,
        validation_json_out,
        validation_csv_out,
        angle_csv_out,
        transport_csv_out,
    ):
        output.parent.mkdir(parents=True, exist_ok=True)

    template = ElasticCdf.from_file(template_path)
    prob_axis, angle_energy_lab = angle_axes(template)
    angle_energy_cm = 0.5 * angle_energy_lab
    sigma_energy_lab = cross_section_axis(template)
    sigma_energy_cm = 0.5 * sigma_energy_lab

    sigma_total = evaluate_model_cm2(MODELS["elastic"], sigma_energy_cm, BOHR_AREA_CM2)
    sigma_mt = evaluate_model_cm2(MODELS["momentum_transfer"], sigma_energy_cm, BOHR_AREA_CM2)

    fits_by_channel = parse_krstic_dcs_markdown(markdown_path)
    elastic_fits = fits_by_channel["Elastic"]
    elastic_fit_by_energy = {float(fit.energy_ev): fit for fit in elastic_fits}
    coeff_json_out.write_text(
        json.dumps(serialize_shared_coefficients(fits_by_channel), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    tabulated_energy_cm, _, tabulated_summary = build_tabulated_inverse_cdf(
        elastic_fits,
        prob_axis,
        theta_floor=THETA_FLOOR_RAD,
        theta_points=THETA_POINTS,
    )
    runtime_angle_grid, runtime_dcs_summary = build_dcs_first_inverse_cdf_by_energy(
        elastic_fits,
        prob_axis,
        angle_energy_cm,
        theta_floor=THETA_FLOOR_RAD,
        theta_points=THETA_POINTS,
    )
    runtime_angle_grid_compat = to_cdf_compatible_angle_order(runtime_angle_grid)
    scattering_angle = runtime_angle_grid_compat.ravel()

    transport = compute_transport_tables_from_sigma_momentum(
        template,
        cross_section=sigma_total,
        sigma_momentum=sigma_mt,
        sigma_energy=sigma_energy_lab,
        xi_points=args.xi_points,
    )
    runtime_angle_transport_ratio = evaluate_angle_transport_ratio(scattering_angle, template)

    apply_blocks(
        template,
        cross_section=sigma_total,
        scattering_angle=scattering_angle,
        reaction_rate=transport["reaction_rate"],
        i_1_0=transport["I_1_0"],
        i_1_1_up=transport["I_1_1_up"],
        i_1_2_up2=transport["I_1_2_up2"],
        sigv_max=float(transport["sigv_max"]),
        angle_min=0.0,
    )
    data_version = (
        "Data Version 3.5 D + D^+ total-elastic cross section from Krstic integral fits "
        "with DCS-first scattering-angle CDF."
    )
    template.write(
        output_path,
        xs_name="dd_00_el_krstic_total",
        description=(
            "Krstic D+ + D total-elastic tables with integral-priority sigma_t/sigma_mt "
            "and DCS-first scattering-angle CDF."
        ),
        data_version=data_version,
    )

    validation_rows: list[dict[str, object]] = []
    for summary in tabulated_summary:
        energy_cm = float(summary["energy_eV_amu"])
        sigma_total_integral = MODELS["elastic"].evaluate_au(energy_cm) * BOHR_AREA_CM2
        sigma_mt_integral = MODELS["momentum_transfer"].evaluate_au(energy_cm) * BOHR_AREA_CM2
        sigma_vi_integral = MODELS["viscosity"].evaluate_au(energy_cm) * BOHR_AREA_CM2
        fit = elastic_fit_by_energy[energy_cm]
        theta_values = np.geomspace(THETA_FLOOR_RAD, np.pi, THETA_POINTS)
        kernel = fit.theta_kernel_au(theta_values)
        sigma_vi_dcs = float(
            np.trapz(np.sin(theta_values) ** 2 * kernel, x=theta_values) * BOHR_AREA_CM2
        )
        integral_ratio = sigma_mt_integral / sigma_total_integral
        validation_rows.append(
            {
                "energy_cm_ev": energy_cm,
                "energy_lab_ev": 2.0 * energy_cm,
                "sigma_t_dcs_cm2": float(summary["sigma_total_cm2"]),
                "sigma_t_integral_cm2": float(sigma_total_integral),
                "sigma_t_relative_error_dcs_vs_integral": float(
                    summary["sigma_total_cm2"] / sigma_total_integral - 1.0
                ),
                "sigma_mt_dcs_cm2": float(summary["sigma_momentum_cm2"]),
                "sigma_mt_integral_cm2": float(sigma_mt_integral),
                "sigma_mt_relative_error_dcs_vs_integral": float(
                    summary["sigma_momentum_cm2"] / sigma_mt_integral - 1.0
                ),
                "sigma_vi_dcs_cm2": sigma_vi_dcs,
                "sigma_vi_integral_cm2": float(sigma_vi_integral),
                "sigma_vi_relative_error_dcs_vs_integral": float(sigma_vi_dcs / sigma_vi_integral - 1.0),
                "transport_ratio_dcs": float(summary["transport_ratio"]),
                "transport_ratio_integral": float(integral_ratio),
                "transport_ratio_relative_error": float(summary["transport_ratio"] / integral_ratio - 1.0),
            }
        )
    write_csv(validation_csv_out, validation_rows)

    runtime_angle_rows: list[dict[str, object]] = []
    for energy_index, energy_lab in enumerate(angle_energy_lab):
        for prob_index, probability in enumerate(prob_axis):
            runtime_angle_rows.append(
                {
                    "energy_lab_ev": float(energy_lab),
                    "energy_cm_ev": float(0.5 * energy_lab),
                    "probability": float(probability),
                    "theta_rad": float(runtime_angle_grid_compat[energy_index, prob_index]),
                }
            )
    write_csv(angle_csv_out, runtime_angle_rows)

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
    write_csv(transport_csv_out, transport_rows)

    runtime_dcs_ratio = np.array(
        [float(summary["transport_ratio"]) for summary in runtime_dcs_summary],
        dtype=float,
    )
    runtime_dcs_cdf_ratio = np.array(
        [float(summary["transport_ratio_from_inverse_cdf"]) for summary in runtime_dcs_summary],
        dtype=float,
    )
    runtime_integral_ratio = np.array(
        [
            MODELS["momentum_transfer"].evaluate_au(float(energy_cm))
            / MODELS["elastic"].evaluate_au(float(energy_cm))
            for energy_cm in angle_energy_cm
        ],
        dtype=float,
    )

    validation = {
        "source_coeff_markdown": str(markdown_path.resolve()),
        "source_coeff_json": coeff_json_out.name,
        "data_version": data_version,
        "cross_section_source": "Krstic integral fit (elastic total)",
        "sigma_momentum_source": "Krstic integral fit (momentum transfer)",
        "angle_source": "Krstic elastic-total direct DCS",
        "angle_interpolation_policy": (
            "DCS-first: interpolate p(theta,E)=2*pi*sin(theta)*d sigma/d Omega "
            "in log(E)-log(p) on the theta grid, then integrate to an inverse CDF."
        ),
        "angle_energy_mapping": "E_cm = 0.5 * E_lab",
        "tabulated_energy_range_cm_ev": [
            float(tabulated_energy_cm[0]),
            float(tabulated_energy_cm[-1]),
        ],
        "runtime_cross_section_range_lab_ev": [
            float(sigma_energy_lab[0]),
            float(sigma_energy_lab[-1]),
        ],
        "runtime_angle_range_lab_ev": [
            float(angle_energy_lab[0]),
            float(angle_energy_lab[-1]),
        ],
        "below_tabulated_energy_policy": "Clamp to the 0.1 eV direct-DCS slice.",
        "theta_floor_rad": THETA_FLOOR_RAD,
        "theta_points": THETA_POINTS,
        "xi_points": args.xi_points,
        "max_abs_sigma_t_relative_error_dcs_vs_integral": float(
            max(abs(row["sigma_t_relative_error_dcs_vs_integral"]) for row in validation_rows)
        ),
        "max_abs_sigma_mt_relative_error_dcs_vs_integral": float(
            max(abs(row["sigma_mt_relative_error_dcs_vs_integral"]) for row in validation_rows)
        ),
        "max_abs_transport_ratio_relative_error": float(
            max(abs(row["transport_ratio_relative_error"]) for row in validation_rows)
        ),
        "max_abs_sigma_vi_relative_error_dcs_vs_integral": float(
            max(abs(row["sigma_vi_relative_error_dcs_vs_integral"]) for row in validation_rows)
        ),
        "max_abs_runtime_angle_vs_dcs_moment_ratio_error": float(
            np.max(np.abs(runtime_angle_transport_ratio - runtime_dcs_ratio))
        ),
        "max_abs_runtime_angle_vs_dcs_inverse_cdf_ratio_error": float(
            np.max(np.abs(runtime_angle_transport_ratio - runtime_dcs_cdf_ratio))
        ),
        "max_abs_runtime_angle_vs_integral_ratio_error": float(
            np.max(np.abs(runtime_angle_transport_ratio - runtime_integral_ratio))
        ),
        "runtime_dcs_rows": runtime_dcs_summary,
        "sigv_max_cm3_s": float(transport["sigv_max"]),
        "validation_rows": validation_rows,
    }
    validation_json_out.write_text(
        json.dumps(validation, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    print(f"Wrote {output_path}")
    print(f"Wrote {validation_json_out}")
    print(f"Wrote {validation_csv_out}")
    print(f"Wrote {angle_csv_out}")
    print(f"Wrote {transport_csv_out}")
    print(f"  sigv_max: {float(transport['sigv_max']):.6e} cm^3/s")


if __name__ == "__main__":
    main()
