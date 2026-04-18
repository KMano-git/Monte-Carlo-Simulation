#!/usr/bin/env python3
"""Build a dd_00_elastic.cdf-compatible CDL file from split CSV sources."""

from __future__ import annotations

import argparse
from pathlib import Path

from cdf_compat import (
    ElasticCdf,
    angle_axes,
    apply_blocks,
    compute_transport_tables,
    load_cross_section_csv,
    load_scalars_json,
    load_scattering_angle_csv,
    load_transport_csv,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--template-cdf", required=True, help="Existing CDF/CDL template to clone")
    parser.add_argument("--output-cdf", required=True, help="Output .cdf/.cdl path")
    parser.add_argument("--cross-section-csv", required=True, help="CSV with cross_section_cm2 column")
    parser.add_argument(
        "--cross-section-column",
        default="cross_section_cm2",
        help="Column name to read from --cross-section-csv",
    )
    parser.add_argument("--scattering-angle-csv", help="CSV with energy/probability/angle_rad rows")
    parser.add_argument(
        "--angle-source-cdf",
        help="Alternate CDF/CDL file to copy scattering_angle from when no angle CSV is provided",
    )
    parser.add_argument(
        "--transport-csv",
        help="Optional CSV with reaction_rate/I_1_x columns; when omitted, transport tables are recomputed",
    )
    parser.add_argument("--scalars-json", help="Optional JSON with sigv_max_cm3_s and angle_min_rad")
    parser.add_argument("--xs-name", help="Replacement xs_name string")
    parser.add_argument("--description", help="Replacement global description attribute")
    parser.add_argument("--data-version", help="Replacement global data_version attribute")
    parser.add_argument(
        "--xi-points",
        type=int,
        default=1000,
        help="Number of xi quadrature points when recomputing transport tables",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    template = ElasticCdf.from_file(args.template_cdf)

    cross_section = load_cross_section_csv(
        args.cross_section_csv,
        value_column=args.cross_section_column,
    )
    if cross_section.size != template.xs_data_inc[0]:
        raise ValueError(
            f"cross_section size mismatch: expected {template.xs_data_inc[0]}, received {cross_section.size}"
        )

    if args.scattering_angle_csv:
        scattering_angle = load_scattering_angle_csv(args.scattering_angle_csv)
    else:
        angle_source_path = args.angle_source_cdf or args.template_cdf
        scattering_angle = ElasticCdf.from_file(angle_source_path).get_block(5)
    if scattering_angle.size != template.xs_data_inc[5]:
        raise ValueError(
            "scattering_angle size mismatch: "
            f"expected {template.xs_data_inc[5]}, received {scattering_angle.size}"
        )

    sigv_max = None
    angle_min = None
    if args.scalars_json:
        scalars = load_scalars_json(args.scalars_json)
        sigv_max = float(scalars.get("sigv_max_cm3_s")) if "sigv_max_cm3_s" in scalars else None
        angle_min = float(scalars.get("angle_min_rad")) if "angle_min_rad" in scalars else None

    if args.transport_csv:
        transport = load_transport_csv(args.transport_csv)
        if sigv_max is None:
            sigv_max = float(transport["reaction_rate"].max())
    else:
        transport = compute_transport_tables(
            template,
            cross_section,
            scattering_angle,
            xi_points=args.xi_points,
        )
        sigv_max = float(transport["sigv_max"])

    if angle_min is None:
        angle_min = float(template.get_block(7)[0])

    apply_blocks(
        template,
        cross_section=cross_section,
        scattering_angle=scattering_angle,
        reaction_rate=transport["reaction_rate"],
        i_1_0=transport["I_1_0"],
        i_1_1_up=transport["I_1_1_up"],
        i_1_2_up2=transport["I_1_2_up2"],
        sigv_max=sigv_max,
        angle_min=angle_min,
    )
    template.write(
        args.output_cdf,
        xs_name=args.xs_name,
        description=args.description,
        data_version=args.data_version,
    )

    _, energy_axis = angle_axes(template)
    print(f"Wrote {args.output_cdf}")
    print(f"  cross_section points : {cross_section.size}")
    print(f"  scattering_angle rows: {scattering_angle.size} ({energy_axis.size} energy slices)")
    print(f"  sigv_max             : {float(template.get_block(6)[0]):.6e} cm^3/s")


if __name__ == "__main__":
    main()
