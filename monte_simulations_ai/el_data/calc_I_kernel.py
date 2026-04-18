#!/usr/bin/env python3
"""Recompute reaction_rate and I_1_x using sigma_mt reconstructed from Krstic DCS."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


def parse_args() -> argparse.Namespace:
    this_dir = Path(__file__).resolve().parent
    out_dir = this_dir / "Krstic"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "input_cdf",
        nargs="?",
        default=str(this_dir / "../test/cdf_change/dd_00_elastic_pure_el_angle.cdf"),
        help="Input CDF/CDL file.",
    )
    parser.add_argument(
        "output_cdf",
        nargs="?",
        default=str(out_dir / "krstic_dd_pure_el_angle_fixed.cdf"),
        help="Output CDF/CDL file.",
    )
    parser.add_argument(
        "--memo",
        default=str(this_dir / "DpD_fit_memo_v2.md"),
        help="Markdown memo holding the manually curated Krstic DCS coefficients.",
    )
    parser.add_argument(
        "--coeff-json",
        default=str(out_dir / "krstic_pure_dcs_coeffs.json"),
        help="Krstic DCS coefficient JSON. It is regenerated automatically when missing.",
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
        help="Log-spaced theta points below the split angle when building sigma_mt from DCS.",
    )
    parser.add_argument(
        "--theta-linear-points",
        type=int,
        default=8192,
        help="Linearly spaced theta points above the split angle when building sigma_mt from DCS.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    from cdf_compat import (
        ElasticCdf,
        angle_axes,
        apply_blocks,
        compute_transport_tables_from_sigma_momentum,
        cross_section_axis,
    )
    from krstic_dcs import (
        build_tabulated_observables,
        interpolate_positive_observable,
        load_or_build_coefficients,
        write_coefficients_json,
    )

    input_path = Path(args.input_cdf).resolve()
    output_path = Path(args.output_cdf).resolve()
    coeff_json_path = Path(args.coeff_json).resolve()
    memo_path = Path(args.memo).resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    coeff_json_path.parent.mkdir(parents=True, exist_ok=True)

    write_coefficients_json(memo_path, coeff_json_path)
    energy_slices = load_or_build_coefficients(
        coeff_json_path=coeff_json_path,
        markdown_path=memo_path,
    )

    cdf = ElasticCdf.from_file(input_path)
    prob_axis, _ = angle_axes(cdf)
    tabulated = build_tabulated_observables(
        energy_slices,
        prob_axis,
        log_points=args.theta_log_points,
        linear_points=args.theta_linear_points,
        negative_policy=args.negative_pure_policy,
    )

    sigma_energy_lab = cross_section_axis(cdf)
    sigma_energy_cm = 0.5 * sigma_energy_lab
    sigma_total = cdf.get_block(0)
    sigma_mt_from_dcs = interpolate_positive_observable(
        tabulated["energy_cm_ev"],
        tabulated["sigma_mt_pure_cm2"],
        sigma_energy_cm,
    )
    transport = compute_transport_tables_from_sigma_momentum(
        cdf,
        cross_section=sigma_total,
        sigma_momentum=sigma_mt_from_dcs,
        sigma_energy=sigma_energy_lab,
        xi_points=args.xi_points,
    )

    apply_blocks(
        cdf,
        cross_section=sigma_total,
        scattering_angle=cdf.get_block(5),
        reaction_rate=transport["reaction_rate"],
        i_1_0=transport["I_1_0"],
        i_1_1_up=transport["I_1_1_up"],
        i_1_2_up2=transport["I_1_2_up2"],
        sigv_max=float(transport["sigv_max"]),
        angle_min=float(cdf.get_block(7)[0]),
    )
    cdf.write(output_path)

    print(f"Wrote {output_path}")
    print(
        "  sigma_mt source               : Krstic DCS (pure elastic, E_cm = 0.5 * E_lab)"
    )
    print(f"  negative pure policy          : {args.negative_pure_policy}")
    print(f"  sigma_mt grid points          : {sigma_mt_from_dcs.size}")
    print(f"  sigv_max                      : {float(transport['sigv_max']):.6e} cm^3/s")


if __name__ == "__main__":
    main()
