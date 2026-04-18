#!/usr/bin/env python3
"""Extract Bachmann-related working data from the repository CDF and memo."""

from __future__ import annotations

import csv
import json
import re
from pathlib import Path


def parse_array(text: str, name: str) -> list[str]:
    pattern = re.compile(rf"\b{re.escape(name)}\s*=\s*(.*?)\s*;", re.DOTALL)
    match = pattern.search(text)
    if not match:
        raise ValueError(f"{name!r} was not found")
    raw = match.group(1).replace("\n", " ")
    return [value.strip() for value in raw.split(",") if value.strip()]


def parse_int_array(text: str, name: str) -> list[int]:
    return [int(value) for value in parse_array(text, name)]


def parse_float_array(text: str, name: str) -> list[float]:
    return [float(value) for value in parse_array(text, name)]


def make_axis(count: int, value_min: float, value_max: float, spacing: str) -> list[float]:
    spacing_clean = spacing.replace('"', "").strip()
    if count == 1:
        return [value_min]
    if spacing_clean == "log" and value_min > 0.0 and value_max > 0.0:
        log_min = math_log(value_min)
        delta = math_log(value_max / value_min) / float(count - 1)
        return [math_exp(log_min + delta * idx) for idx in range(count)]
    step = (value_max - value_min) / float(count - 1)
    return [value_min + step * idx for idx in range(count)]


def math_log(value: float) -> float:
    import math

    return math.log(value)


def math_exp(value: float) -> float:
    import math

    return math.exp(value)


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"No rows were provided for {path}")
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def extract_memo_excerpt(memo_text: str) -> str:
    start_tag = "## Bachmann PDF page 52"
    end_tag = "## Bachmann PDF page 55"
    start = memo_text.find(start_tag)
    end = memo_text.find(end_tag)
    if start < 0:
        raise ValueError("Bachmann memo excerpt start marker was not found")
    if end < 0:
        raise ValueError("Bachmann memo excerpt end marker was not found")
    excerpt = memo_text[start:end].rstrip() + "\n"
    return excerpt


def main() -> None:
    out_dir = Path(__file__).resolve().parent
    el_data_dir = out_dir.parent

    cdf_path = el_data_dir / "dd_00_elastic.cdf"
    memo_path = el_data_dir / "DpD_fit_memo.md"

    cdf_text = cdf_path.read_text(encoding="utf-8")
    memo_text = memo_path.read_text(encoding="utf-8")

    xs_data_base = parse_int_array(cdf_text, "xs_data_base")
    xs_data_inc = parse_int_array(cdf_text, "xs_data_inc")
    xs_data_tab = parse_float_array(cdf_text, "xs_data_tab")

    xs_tab_index_flat = parse_int_array(cdf_text, "xs_tab_index")
    xs_tab_index = [xs_tab_index_flat[idx * 3 : (idx + 1) * 3] for idx in range(20)]

    xs_min_flat = parse_float_array(cdf_text, "xs_min")
    xs_min = [xs_min_flat[idx * 3 : (idx + 1) * 3] for idx in range(20)]

    xs_max_flat = parse_float_array(cdf_text, "xs_max")
    xs_max = [xs_max_flat[idx * 3 : (idx + 1) * 3] for idx in range(20)]

    xs_spacing_flat = parse_array(cdf_text, "xs_spacing")
    xs_spacing = [xs_spacing_flat[idx * 4 : (idx + 1) * 4] for idx in range(20)]

    def get_data(block_index: int) -> list[float]:
        base = xs_data_base[block_index]
        size = xs_data_inc[block_index]
        return xs_data_tab[base : base + size]

    cross_section = get_data(0)
    reaction_rate = get_data(1)
    i_1_0 = get_data(2)
    i_1_1_up = get_data(3)
    i_1_2_up2 = get_data(4)
    scattering_angle = get_data(5)
    sigv_max = get_data(6)[0]
    angle_min = get_data(7)[0]

    energy_sigma = make_axis(
        xs_tab_index[0][0], xs_min[0][0], xs_max[0][0], xs_spacing[0][0]
    )
    energy_transport = make_axis(
        xs_tab_index[1][0], xs_min[1][0], xs_max[1][0], xs_spacing[1][0]
    )
    temp_transport = make_axis(
        xs_tab_index[1][1], xs_min[1][1], xs_max[1][1], xs_spacing[1][1]
    )
    probability_axis = make_axis(
        xs_tab_index[5][0], xs_min[5][0], xs_max[5][0], xs_spacing[5][0]
    )
    energy_angle = make_axis(
        xs_tab_index[5][1], xs_min[5][1], xs_max[5][1], xs_spacing[5][1]
    )

    cross_rows = [
        {
            "energy_eV_amu": energy_sigma[idx],
            "cross_section_cm2": cross_section[idx],
        }
        for idx in range(len(cross_section))
    ]

    transport_rows: list[dict[str, object]] = []
    n_energy_transport = xs_tab_index[1][0]
    n_temp_transport = xs_tab_index[1][1]
    for temp_idx in range(n_temp_transport):
        for energy_idx in range(n_energy_transport):
            flat_index = temp_idx * n_energy_transport + energy_idx
            transport_rows.append(
                {
                    "energy_eV_amu": energy_transport[energy_idx],
                    "temperature_eV_amu": temp_transport[temp_idx],
                    "reaction_rate_cm3_s": reaction_rate[flat_index],
                    "I_1_0_cm3_s": i_1_0[flat_index],
                    "I_1_1_up_cm4_s2": i_1_1_up[flat_index],
                    "I_1_2_up2_cm5_s3": i_1_2_up2[flat_index],
                }
            )

    angle_rows: list[dict[str, object]] = []
    n_probability = xs_tab_index[5][0]
    n_energy_angle = xs_tab_index[5][1]
    for energy_idx in range(n_energy_angle):
        for prob_idx in range(n_probability):
            flat_index = energy_idx * n_probability + prob_idx
            angle_rows.append(
                {
                    "energy_eV_amu": energy_angle[energy_idx],
                    "probability": probability_axis[prob_idx],
                    "angle_rad": scattering_angle[flat_index],
                }
            )

    scalars = {
        "source_cdf": str(cdf_path.relative_to(el_data_dir.parent)),
        "sigv_max_cm3_s": sigv_max,
        "angle_min_rad": angle_min,
        "cross_section_points": len(cross_rows),
        "transport_rows": len(transport_rows),
        "angle_rows": len(angle_rows),
    }

    write_csv(out_dir / "bachmann_dd_cross_section_from_cdf.csv", cross_rows)
    write_csv(out_dir / "bachmann_dd_transport_tables_from_cdf.csv", transport_rows)
    write_csv(out_dir / "bachmann_dd_scattering_angle_from_cdf.csv", angle_rows)
    (out_dir / "bachmann_dd_scalars_from_cdf.json").write_text(
        json.dumps(scalars, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    (out_dir / "bachmann_raw_coefficients_excerpt.md").write_text(
        extract_memo_excerpt(memo_text),
        encoding="utf-8",
    )

    print("Wrote bachmann_dd_cross_section_from_cdf.csv")
    print("Wrote bachmann_dd_transport_tables_from_cdf.csv")
    print("Wrote bachmann_dd_scattering_angle_from_cdf.csv")
    print("Wrote bachmann_dd_scalars_from_cdf.json")
    print("Wrote bachmann_raw_coefficients_excerpt.md")


if __name__ == "__main__":
    main()
