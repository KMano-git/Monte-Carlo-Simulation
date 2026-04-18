#!/usr/bin/env python3
"""Generate structured D+ + D integral cross-section datasets from Krstic fits."""

from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path


BOHR_AREA_CM2 = 2.80028e-17
FIT_RANGE_MIN_EV = 0.1
FIT_RANGE_MAX_EV = 100.0


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


@dataclass(frozen=True)
class LogLogInterpolator:
    label: str
    points: tuple[tuple[float, float], ...]

    def evaluate_au(self, energy_ev: float) -> float:
        if energy_ev <= 0.0:
            raise ValueError(f"Energy must be positive for {self.label}")

        log_energy = math.log(energy_ev)
        log_points = [(math.log(energy), math.log(value)) for energy, value in self.points]

        if log_energy <= log_points[0][0]:
            x0, y0 = log_points[0]
            x1, y1 = log_points[1]
        elif log_energy >= log_points[-1][0]:
            x0, y0 = log_points[-2]
            x1, y1 = log_points[-1]
        else:
            x0 = y0 = x1 = y1 = 0.0
            for (left_x, left_y), (right_x, right_y) in zip(log_points[:-1], log_points[1:]):
                if left_x <= log_energy <= right_x:
                    x0, y0 = left_x, left_y
                    x1, y1 = right_x, right_y
                    break

        if abs(x1 - x0) < 1.0e-30:
            return math.exp(y0)

        weight = (log_energy - x0) / (x1 - x0)
        return math.exp((1.0 - weight) * y0 + weight * y1)


ANALYTIC_MODELS = {
    "elastic": RationalLogFit(
        label="elastic",
        a=(0.637254e03, -0.125310e03, 0.167078e02, -0.990209e00),
    ),
    "momentum_transfer": RationalLogFit(
        label="momentum_transfer",
        a=(0.347935e03, 0.419075e03, 0.117426e03, -0.324821e01, -0.105819e01),
        b=(0.131735e01, 0.482389e00, 0.401043e-01),
    ),
}


REFERENCE_POINTS = [
    {
        "energy_eV_cm": 0.1000,
        "elastic_au": 9.773961e02,
        "momentum_transfer_au": 4.120358e02,
        "viscosity_au": 1.813288e02,
        "spin_exchange_au": 2.213106e02,
    },
    {
        "energy_eV_cm": 0.1995,
        "elastic_au": 9.096440e02,
        "momentum_transfer_au": 4.008971e02,
        "viscosity_au": 1.358733e02,
        "spin_exchange_au": 2.056762e02,
    },
    {
        "energy_eV_cm": 0.5012,
        "elastic_au": 7.024936e02,
        "momentum_transfer_au": 3.630347e02,
        "viscosity_au": 1.079634e02,
        "spin_exchange_au": 1.802378e02,
    },
    {
        "energy_eV_cm": 1.0000,
        "elastic_au": 6.461367e02,
        "momentum_transfer_au": 3.487627e02,
        "viscosity_au": 9.545577e01,
        "spin_exchange_au": 1.752253e02,
    },
    {
        "energy_eV_cm": 1.9950,
        "elastic_au": 5.621003e02,
        "momentum_transfer_au": 3.241112e02,
        "viscosity_au": 6.752089e01,
        "spin_exchange_au": 1.617602e02,
    },
    {
        "energy_eV_cm": 5.0120,
        "elastic_au": 4.986651e02,
        "momentum_transfer_au": 2.897988e02,
        "viscosity_au": 2.707555e01,
        "spin_exchange_au": 1.448380e02,
    },
    {
        "energy_eV_cm": 10.0000,
        "elastic_au": 4.488888e02,
        "momentum_transfer_au": 2.654433e02,
        "viscosity_au": 1.410217e01,
        "spin_exchange_au": 1.326478e02,
    },
    {
        "energy_eV_cm": 19.9500,
        "elastic_au": 3.923312e02,
        "momentum_transfer_au": 2.419854e02,
        "viscosity_au": 7.248538e00,
        "spin_exchange_au": 1.209559e02,
    },
    {
        "energy_eV_cm": 50.1200,
        "elastic_au": 3.492420e02,
        "momentum_transfer_au": 2.121037e02,
        "viscosity_au": 2.436365e00,
        "spin_exchange_au": 1.060618e02,
    },
    {
        "energy_eV_cm": 100.0000,
        "elastic_au": 3.011882e02,
        "momentum_transfer_au": 1.907672e02,
        "viscosity_au": 8.959529e-01,
        "spin_exchange_au": 9.538493e01,
    },
]


def build_reference_lookup(key: str) -> tuple[tuple[float, float], ...]:
    return tuple((float(point["energy_eV_cm"]), float(point[f"{key}_au"])) for point in REFERENCE_POINTS)


MODELS = {
    **ANALYTIC_MODELS,
    "viscosity": LogLogInterpolator(
        label="viscosity",
        points=build_reference_lookup("viscosity"),
    ),
    "spin_exchange": LogLogInterpolator(
        label="spin_exchange",
        points=build_reference_lookup("spin_exchange"),
    ),
}


def to_serializable_coeffs() -> dict[str, object]:
    fits_json: dict[str, object] = {}
    for key, model in MODELS.items():
        if isinstance(model, RationalLogFit):
            fits_json[key] = {
                "method": "analytic_rational_log_fit",
                "a_coeffs": list(model.a),
                "b_coeffs": list(model.b),
                "formula": "sigma(E) = (sum_i a_i (ln E)^i) / (1 + sum_j b_j (ln E)^(j+1))",
                "units_input_energy": "eV (center of mass)",
                "units_output_au": "a0^2",
                "units_output_cm2": "cm^2",
            }
        else:
            fits_json[key] = {
                "method": "loglog_interpolation_from_page104_reference_points",
                "reason": (
                    "The OCR-transcribed analytic coefficients are incomplete or unreliable "
                    "for this channel in the current memo snapshot, so the structured working "
                    "dataset is generated from the tabulated page-104 reference values."
                ),
                "reference_points": [
                    {"energy_eV_cm": energy, "sigma_au": sigma}
                    for energy, sigma in model.points
                ],
                "units_input_energy": "eV (center of mass)",
                "units_output_au": "a0^2",
                "units_output_cm2": "cm^2",
            }
    return {
        "source": "monte_simulations_ai/el_data/DpD_fit_memo.md",
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


if __name__ == "__main__":
    main()
