#!/usr/bin/env python3
"""Build a direct-DCS Krstic composite dataset in dd_00_elastic.cdf layout."""

from __future__ import annotations

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


def main() -> None:
    this_dir = Path(__file__).resolve().parent
    el_data_dir = this_dir.parent
    repo_dir = el_data_dir.parent
    sys.path.insert(0, str(el_data_dir))
    sys.path.insert(0, str(this_dir))

    from cdf_compat import (
        ElasticCdf,
        angle_axes,
        apply_blocks,
        compute_transport_tables,
        evaluate_angle_transport_ratio,
        load_cross_section_csv,
        transport_axes,
    )
    from generate_krstic_integral_data import BOHR_AREA_CM2, MODELS
    from krstic_dcs import (
        build_tabulated_inverse_cdf,
        interpolate_inverse_cdf_by_energy,
        parse_krstic_dcs_markdown,
        serialize_fits,
    )

    template_path = el_data_dir / "dd_00_elastic.cdf"
    markdown_path = el_data_dir / "DpD_fit_memo_v2.md"
    output_path = this_dir / "krstic_dd_composite_compat.cdf"
    coeff_json_out = this_dir / "krstic_dd_dcs_coeffs.json"
    dcs_checks_csv_out = this_dir / "krstic_dd_dcs_integral_checks.csv"
    tabulated_angle_csv_out = this_dir / "krstic_dd_scattering_angle_tabulated_31.csv"
    runtime_angle_csv_out = this_dir / "krstic_dd_scattering_angle_compat.csv"
    transport_csv_out = this_dir / "krstic_dd_transport_compat.csv"
    validation_json_out = this_dir / "krstic_dd_angle_validation.json"
    cross_section_csv = this_dir / "krstic_dd_integral_cdf_grid_101.csv"

    template = ElasticCdf.from_file(template_path)
    cross_section = load_cross_section_csv(cross_section_csv, value_column="elastic_cm2")
    prob_axis, runtime_angle_energy = angle_axes(template)

    fits_by_channel = parse_krstic_dcs_markdown(markdown_path)
    elastic_fits = fits_by_channel["Elastic"]
    coeff_json_out.write_text(
        json.dumps(
            serialize_fits(
                fits_by_channel,
                source_markdown=str(markdown_path.relative_to(repo_dir)),
            ),
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    tabulated_energy, tabulated_angle_table, tabulated_summary = build_tabulated_inverse_cdf(
        elastic_fits,
        prob_axis,
        theta_floor=THETA_FLOOR_RAD,
        theta_points=THETA_POINTS,
    )
    runtime_angle_grid = interpolate_inverse_cdf_by_energy(
        tabulated_energy,
        tabulated_angle_table,
        runtime_angle_energy,
    )
    runtime_angle_grid_compat = to_cdf_compatible_angle_order(runtime_angle_grid)
    scattering_angle = runtime_angle_grid_compat.ravel()

    transport = compute_transport_tables(template, cross_section, scattering_angle)
    runtime_transport_ratio = evaluate_angle_transport_ratio(scattering_angle, template)

    apply_blocks(
        template,
        cross_section=cross_section,
        scattering_angle=scattering_angle,
        reaction_rate=transport["reaction_rate"],
        i_1_0=transport["I_1_0"],
        i_1_1_up=transport["I_1_1_up"],
        i_1_2_up2=transport["I_1_2_up2"],
        sigv_max=float(transport["sigv_max"]),
        angle_min=0.0,
    )
    template.write(
        output_path,
        xs_name="dd_00_el_krstic_comp",
        description="Krstic composite D + D^+ elastic tables with direct differential-fit angles.",
        data_version=(
            "Data Version 3.2 D + D^+ Krstic composite cross section with direct "
            "differential-fit angle CDF."
        ),
    )

    tabulated_angle_rows: list[dict[str, object]] = []
    for energy_index, energy_value in enumerate(tabulated_energy):
        for prob_index, prob_value in enumerate(prob_axis):
            tabulated_angle_rows.append(
                {
                    "energy_eV_amu": float(energy_value),
                    "probability": float(prob_value),
                    "angle_rad": float(tabulated_angle_table[energy_index, prob_index]),
                    "source": "direct_dcs",
                }
            )
    write_csv(tabulated_angle_csv_out, tabulated_angle_rows)

    runtime_angle_rows: list[dict[str, object]] = []
    for energy_index, energy_value in enumerate(runtime_angle_energy):
        for prob_index, prob_value in enumerate(prob_axis):
            runtime_angle_rows.append(
                {
                    "energy_eV_amu": float(energy_value),
                    "probability": float(prob_value),
                    "angle_rad": float(runtime_angle_grid_compat[energy_index, prob_index]),
                }
            )
    write_csv(runtime_angle_csv_out, runtime_angle_rows)

    transport_rows_out: list[dict[str, object]] = []
    energy_transport, temp_transport = transport_axes(template)
    rr_table = transport["reaction_rate"].reshape(temp_transport.size, energy_transport.size)
    i10_table = transport["I_1_0"].reshape(temp_transport.size, energy_transport.size)
    i11_table = transport["I_1_1_up"].reshape(temp_transport.size, energy_transport.size)
    i12_table = transport["I_1_2_up2"].reshape(temp_transport.size, energy_transport.size)
    for temp_index, temp_value in enumerate(temp_transport):
        for energy_index, energy_value in enumerate(energy_transport):
            transport_rows_out.append(
                {
                    "energy_eV_amu": float(energy_value),
                    "temperature_eV_amu": float(temp_value),
                    "reaction_rate_cm3_s": float(rr_table[temp_index, energy_index]),
                    "I_1_0_cm3_s": float(i10_table[temp_index, energy_index]),
                    "I_1_1_up_cm4_s2": float(i11_table[temp_index, energy_index]),
                    "I_1_2_up2_cm5_s3": float(i12_table[temp_index, energy_index]),
                }
            )
    write_csv(transport_csv_out, transport_rows_out)

    dcs_check_rows: list[dict[str, object]] = []
    tabulated_integral_ratio = np.zeros(tabulated_energy.size, dtype=float)
    tabulated_direct_ratio = np.zeros(tabulated_energy.size, dtype=float)
    for index, summary in enumerate(tabulated_summary):
        energy_value = float(summary["energy_eV_amu"])
        sigma_total_integral = MODELS["elastic"].evaluate_au(energy_value) * BOHR_AREA_CM2
        sigma_momentum_integral = (
            MODELS["momentum_transfer"].evaluate_au(energy_value) * BOHR_AREA_CM2
        )
        integral_ratio = sigma_momentum_integral / sigma_total_integral
        tabulated_integral_ratio[index] = integral_ratio
        tabulated_direct_ratio[index] = float(summary["transport_ratio"])
        dcs_check_rows.append(
            {
                "energy_eV_amu": energy_value,
                "sigma_total_dcs_cm2": float(summary["sigma_total_cm2"]),
                "sigma_total_integral_cm2": float(sigma_total_integral),
                "sigma_total_relative_error": float(summary["sigma_total_cm2"] / sigma_total_integral - 1.0),
                "sigma_momentum_dcs_cm2": float(summary["sigma_momentum_cm2"]),
                "sigma_momentum_integral_cm2": float(sigma_momentum_integral),
                "sigma_momentum_relative_error": float(
                    summary["sigma_momentum_cm2"] / sigma_momentum_integral - 1.0
                ),
                "sigma_m_over_sigma_t_dcs": float(summary["transport_ratio"]),
                "sigma_m_over_sigma_t_integral": float(integral_ratio),
                "transport_ratio_relative_error": float(summary["transport_ratio"] / integral_ratio - 1.0),
            }
        )
    write_csv(dcs_checks_csv_out, dcs_check_rows)

    interpolated_direct_ratio = np.interp(
        np.log(np.clip(runtime_angle_energy, tabulated_energy[0], tabulated_energy[-1])),
        np.log(tabulated_energy),
        tabulated_direct_ratio,
        left=float(tabulated_direct_ratio[0]),
        right=float(tabulated_direct_ratio[-1]),
    )
    runtime_integral_ratio = np.array(
        [
            MODELS["momentum_transfer"].evaluate_au(float(energy_value))
            / MODELS["elastic"].evaluate_au(float(energy_value))
            for energy_value in runtime_angle_energy
        ],
        dtype=float,
    )

    validation = {
        "source_coeff_markdown": str(markdown_path.relative_to(repo_dir)),
        "source_coeff_json": coeff_json_out.name,
        "source_cross_section_csv": cross_section_csv.name,
        "source_angle_shape": "direct_krstic_dcs_logE_interp_inverse_cdf",
        "tabulated_energy_range_eV_amu": [float(tabulated_energy[0]), float(tabulated_energy[-1])],
        "runtime_energy_range_eV_amu": [
            float(runtime_angle_energy[0]),
            float(runtime_angle_energy[-1]),
        ],
        "below_tabulated_energy_policy": "Clamp to the 0.1 eV direct-DCS slice.",
        "theta_floor_rad": THETA_FLOOR_RAD,
        "theta_points": THETA_POINTS,
        "max_abs_tabulated_sigma_total_relative_error": float(
            np.max(np.abs([row["sigma_total_relative_error"] for row in dcs_check_rows]))
        ),
        "max_abs_tabulated_transport_ratio_relative_error": float(
            np.max(np.abs([row["transport_ratio_relative_error"] for row in dcs_check_rows]))
        ),
        "max_abs_runtime_vs_interpolated_direct_ratio_error": float(
            np.max(np.abs(runtime_transport_ratio - interpolated_direct_ratio))
        ),
        "max_abs_runtime_vs_integral_ratio_error": float(
            np.max(np.abs(runtime_transport_ratio - runtime_integral_ratio))
        ),
        "sigv_max_cm3_s": float(transport["sigv_max"]),
        "tabulated_per_energy": dcs_check_rows,
        "runtime_per_energy": [
            {
                "energy_eV_amu": float(runtime_angle_energy[index]),
                "angle_table_transport_ratio": float(runtime_transport_ratio[index]),
                "interpolated_direct_transport_ratio": float(interpolated_direct_ratio[index]),
                "integral_transport_ratio": float(runtime_integral_ratio[index]),
            }
            for index in range(runtime_angle_energy.size)
        ],
    }
    validation_json_out.write_text(
        json.dumps(validation, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    print(f"Wrote {coeff_json_out}")
    print(f"Wrote {dcs_checks_csv_out}")
    print(f"Wrote {tabulated_angle_csv_out}")
    print(f"Wrote {runtime_angle_csv_out}")
    print(f"Wrote {transport_csv_out}")
    print(f"Wrote {validation_json_out}")
    print(f"Wrote {output_path}")
    print(f"  sigv_max                              : {float(transport['sigv_max']):.6e} cm^3/s")
    print(
        "  max |runtime ratio - interpolated direct| : "
        f"{float(np.max(np.abs(runtime_transport_ratio - interpolated_direct_ratio))):.6e}"
    )
    print(
        "  max |runtime ratio - integral ratio|      : "
        f"{float(np.max(np.abs(runtime_transport_ratio - runtime_integral_ratio))):.6e}"
    )


if __name__ == "__main__":
    main()
