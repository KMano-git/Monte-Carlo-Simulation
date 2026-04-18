#!/usr/bin/env python3
"""Generate structured D+ + D integral cross-section datasets from Krstic fits."""

from __future__ import annotations

import ast
import csv
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path


BOHR_AREA_CM2 = 2.80028e-17
FIT_RANGE_MIN_EV = 0.1
FIT_RANGE_MAX_EV = 100.0
EL_DATA_DIR = Path(__file__).resolve().parent.parent
DEFAULT_MEMO_PATH = EL_DATA_DIR / "DpD_fit_memo_v2.md"
REFERENCE_TABLE_SECTION_MARKER = "## 2. Krstić: D⁺ + D 積分断面積表（OCR補正反映済み）"
MANUAL_INTEGRAL_SECTION_MARKER = "## 12. Cross Sectionのfitting parametersの手動入力欄"
CHANNEL_NAME_MAP = {
    "Elastic": "elastic",
    "Momentum Transfer": "momentum_transfer",
    "Viscosity": "viscosity",
    "Spin Exchange": "spin_exchange",
}


def build_log_grid(count: int, energy_min: float, energy_max: float) -> list[float]:
    log_min = math.log(energy_min)
    delta = math.log(energy_max / energy_min) / float(count - 1)
    return [math.exp(log_min + delta * idx) for idx in range(count)]


CDF_GRID_101 = build_log_grid(101, 9.9318169312296e-4, 9.9318169312296e1)
ANGLE_GRID_51 = build_log_grid(51, 9.9318169312296e-4, 9.9318169312296e1)


@dataclass(frozen=True)
class RationalLogFit:
    label: str
    a: tuple[float, ...]
    b: tuple[float, ...] = ()

    def evaluate_au(self, energy_ev: float) -> float:
        x_val = math.log(energy_ev)
        numerator = sum(coeff * (x_val ** power) for power, coeff in enumerate(self.a))
        denominator = 1.0 + sum(
            coeff * (x_val ** (power + 1)) for power, coeff in enumerate(self.b)
        )
        if abs(denominator) < 1.0e-30:
            raise ZeroDivisionError(
                f"Denominator nearly vanished for {self.label} at E={energy_ev:.6e} eV"
            )
        sigma_au = numerator / denominator
        if sigma_au <= 0.0:
            raise ValueError(
                f"Non-positive cross-section from analytic fit for {self.label} at "
                f"E={energy_ev:.6e} eV"
            )
        return sigma_au


def _parse_float_list(raw: str) -> tuple[float, ...]:
    return tuple(float(value) for value in ast.literal_eval(raw.replace("E", "e")))


def parse_reference_points(markdown_path: str | Path) -> list[dict[str, float]]:
    text = Path(markdown_path).read_text(encoding="utf-8")
    if REFERENCE_TABLE_SECTION_MARKER not in text:
        raise ValueError(f"{REFERENCE_TABLE_SECTION_MARKER!r} was not found in {markdown_path}")

    body = text.split(REFERENCE_TABLE_SECTION_MARKER, maxsplit=1)[1]
    rows: list[dict[str, float]] = []
    table_started = False
    for line in body.splitlines():
        stripped = line.strip()
        if stripped.startswith("|"):
            table_started = True
            cells = [cell.strip() for cell in stripped.strip("|").split("|")]
            if not cells or cells[0].startswith("---") or cells[0] == "Energy (CM) [eV]":
                continue
            if len(cells) != 5:
                raise ValueError(f"Unexpected reference-table row in {markdown_path}: {line!r}")
            rows.append(
                {
                    "energy_eV_cm": float(cells[0]),
                    "elastic_au": float(cells[1]),
                    "momentum_transfer_au": float(cells[2]),
                    "viscosity_au": float(cells[3]),
                    "spin_exchange_au": float(cells[4]),
                }
            )
            continue
        if table_started:
            break

    if not rows:
        raise ValueError(f"No reference-point rows were parsed from {markdown_path}")
    return rows


def parse_integral_fit_models(markdown_path: str | Path) -> dict[str, RationalLogFit]:
    text = Path(markdown_path).read_text(encoding="utf-8")
    if MANUAL_INTEGRAL_SECTION_MARKER not in text:
        raise ValueError(f"{MANUAL_INTEGRAL_SECTION_MARKER!r} was not found in {markdown_path}")

    body = text.split(MANUAL_INTEGRAL_SECTION_MARKER, maxsplit=1)[1]
    parsed: dict[str, dict[str, tuple[float, ...]]] = {}
    current_channel: str | None = None

    for line in body.splitlines():
        channel_match = re.match(r"^###\s+(Elastic|Momentum Transfer|Viscosity|Spin Exchange)\s*$", line)
        if channel_match:
            current_channel = CHANNEL_NAME_MAP[channel_match.group(1)]
            parsed[current_channel] = {}
            continue

        value_match = re.match(r"^\s*-\s+([A-Za-z_]+):\s*(.*)$", line)
        if value_match and current_channel is not None:
            key = value_match.group(1)
            raw_value = value_match.group(2).strip()
            if key.endswith("_values"):
                parsed[current_channel][key] = _parse_float_list(raw_value)

    models: dict[str, RationalLogFit] = {}
    for channel_key in CHANNEL_NAME_MAP.values():
        channel_data = parsed.get(channel_key)
        if channel_data is None:
            raise ValueError(f"Missing manual integral-fit entry for {channel_key}")
        if "a_values" not in channel_data:
            raise ValueError(f"Missing a_values for {channel_key}")
        models[channel_key] = RationalLogFit(
            label=channel_key,
            a=channel_data["a_values"],
            b=channel_data.get("b_values", ()),
        )
    return models


REFERENCE_POINTS = parse_reference_points(DEFAULT_MEMO_PATH)
MODELS = parse_integral_fit_models(DEFAULT_MEMO_PATH)


def to_serializable_coeffs() -> dict[str, object]:
    fits_json: dict[str, object] = {}
    for key, model in MODELS.items():
        fits_json[key] = {
            "method": "analytic_rational_log_fit_from_manual_entry",
            "a_coeffs": list(model.a),
            "b_coeffs": list(model.b),
            "formula": "sigma(E) = (sum_i a_i (ln E)^i) / (1 + sum_j b_j (ln E)^(j+1))",
            "source_markdown": "monte_simulations_ai/el_data/DpD_fit_memo_v2.md",
            "source_section": MANUAL_INTEGRAL_SECTION_MARKER,
            "units_input_energy": "eV (center of mass)",
            "units_output_au": "a0^2",
            "units_output_cm2": "cm^2",
        }
    return {
        "source": "monte_simulations_ai/el_data/DpD_fit_memo_v2.md",
        "source_reference_section": REFERENCE_TABLE_SECTION_MARKER,
        "source_fit_section": MANUAL_INTEGRAL_SECTION_MARKER,
        "fit_range_eV_cm": [FIT_RANGE_MIN_EV, FIT_RANGE_MAX_EV],
        "bohr_area_cm2": BOHR_AREA_CM2,
        "fits": fits_json,
    }


def build_reference_rows() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for point in REFERENCE_POINTS:
        energy = float(point["energy_eV_cm"])
        row: dict[str, object] = {
            "energy_eV_cm": energy,
            "in_tabulated_fit_range": True,
        }
        for key, model in MODELS.items():
            tabulated_au = float(point[f"{key}_au"])
            fitted_au = model.evaluate_au(energy)
            row[f"{key}_tabulated_au"] = tabulated_au
            row[f"{key}_tabulated_cm2"] = tabulated_au * BOHR_AREA_CM2
            row[f"{key}_fit_au"] = fitted_au
            row[f"{key}_fit_cm2"] = fitted_au * BOHR_AREA_CM2
            row[f"{key}_relative_error"] = (fitted_au / tabulated_au) - 1.0
        rows.append(row)
    return rows


def build_grid_rows(energies: list[float]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for energy in energies:
        row: dict[str, object] = {
            "energy_eV_amu": energy,
            "energy_eV_cm": energy,
            "in_tabulated_fit_range": FIT_RANGE_MIN_EV <= energy <= FIT_RANGE_MAX_EV,
        }
        for key, model in MODELS.items():
            sigma_au = model.evaluate_au(energy)
            row[f"{key}_au"] = sigma_au
            row[f"{key}_cm2"] = sigma_au * BOHR_AREA_CM2
        rows.append(row)
    return rows


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"No rows were provided for {path}")
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    out_dir = Path(__file__).resolve().parent

    coeffs_path = out_dir / "krstic_dd_integral_coeffs.json"
    reference_csv = out_dir / "krstic_dd_integral_reference_points.csv"
    cdf_grid_csv = out_dir / "krstic_dd_integral_cdf_grid_101.csv"
    angle_grid_csv = out_dir / "krstic_dd_integral_angle_grid_51.csv"

    coeffs_path.write_text(
        json.dumps(to_serializable_coeffs(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_csv(reference_csv, build_reference_rows())
    write_csv(cdf_grid_csv, build_grid_rows(CDF_GRID_101))
    write_csv(angle_grid_csv, build_grid_rows(ANGLE_GRID_51))

    print(f"Wrote {coeffs_path.name}")
    print(f"Wrote {reference_csv.name}")
    print(f"Wrote {cdf_grid_csv.name}")
    print(f"Wrote {angle_grid_csv.name}")
    print(f"  memo source: {DEFAULT_MEMO_PATH}")


if __name__ == "__main__":
    main()
