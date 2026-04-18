#!/usr/bin/env python3
"""Rebuild the current Bachmann-based elastic CDF from extracted split tables."""

from __future__ import annotations

import sys
from pathlib import Path


def main() -> None:
    this_dir = Path(__file__).resolve().parent
    el_data_dir = this_dir.parent
    sys.path.insert(0, str(el_data_dir))

    from cdf_compat import (
        ElasticCdf,
        apply_blocks,
        load_cross_section_csv,
        load_scalars_json,
        load_scattering_angle_csv,
        load_transport_csv,
    )

    template_path = el_data_dir / "dd_00_elastic.cdf"
    output_path = this_dir / "bachmann_dd_from_split_tables_compat.cdf"

    template = ElasticCdf.from_file(template_path)
    cross_section = load_cross_section_csv(this_dir / "bachmann_dd_cross_section_from_cdf.csv")
    scattering_angle = load_scattering_angle_csv(this_dir / "bachmann_dd_scattering_angle_from_cdf.csv")
    transport = load_transport_csv(this_dir / "bachmann_dd_transport_tables_from_cdf.csv")
    scalars = load_scalars_json(this_dir / "bachmann_dd_scalars_from_cdf.json")

    apply_blocks(
        template,
        cross_section=cross_section,
        scattering_angle=scattering_angle,
        reaction_rate=transport["reaction_rate"],
        i_1_0=transport["I_1_0"],
        i_1_1_up=transport["I_1_1_up"],
        i_1_2_up2=transport["I_1_2_up2"],
        sigv_max=float(scalars["sigv_max_cm3_s"]),
        angle_min=float(scalars["angle_min_rad"]),
    )
    template.write(
        output_path,
        xs_name="dd_00_el_bachmann_src",
        description="Bachmann-based D + D^+ elastic tables regenerated from extracted split tables.",
        data_version="Data Version 3.1 D + D^+ Bachmann split-table reconstruction.",
    )

    print(f"Wrote {output_path}")
    print("  source tables        : Bachmann/*.csv + bachmann_dd_scalars_from_cdf.json")
    print(f"  sigv_max             : {float(scalars['sigv_max_cm3_s']):.6e} cm^3/s")


if __name__ == "__main__":
    main()
