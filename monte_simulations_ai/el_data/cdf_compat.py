#!/usr/bin/env python3
"""Utilities for reading, updating, and regenerating elastic CDL/CDF tables."""

from __future__ import annotations

import csv
import json
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy import integrate
from scipy.interpolate import interp1d


C0_CM_PER_S = 1.3891494e6
NUM_DEP_VAR = 20
RANK_IND = 3
RANK_IND0 = 4


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


def make_axis(count: int, value_min: float, value_max: float, spacing: str) -> np.ndarray:
    spacing_clean = spacing.replace('"', "").strip()
    if count == 1:
        return np.array([value_min], dtype=float)
    if spacing_clean == "log" and value_min > 0.0 and value_max > 0.0:
        return np.logspace(np.log10(value_min), np.log10(value_max), count)
    return np.linspace(value_min, value_max, count)


def format_float_block(values: np.ndarray, per_line: int = 4) -> str:
    chunks = [f"{float(value):.14e}" for value in np.asarray(values, dtype=float).ravel()]
    lines = []
    for index in range(0, len(chunks), per_line):
        lines.append("    " + ", ".join(chunks[index : index + per_line]))
    return ",\n".join(lines)


def replace_array_block(text: str, name: str, formatted_values: str) -> str:
    pattern = re.compile(rf"(\b{re.escape(name)}\s*=\s*)(.*?)(\s*;)", re.DOTALL)
    match = pattern.search(text)
    if not match:
        raise ValueError(f"{name!r} assignment was not found")
    return pattern.sub(rf"\1\n{formatted_values}\n  \3", text, count=1)


def replace_char_scalar(text: str, name: str, value: str) -> str:
    pattern = re.compile(rf"(\b{re.escape(name)}\s*=\s*)\".*?\"(\s*;)", re.DOTALL)
    if not pattern.search(text):
        raise ValueError(f"{name!r} string assignment was not found")
    return pattern.sub(rf'\1"{value}"\2', text, count=1)


def replace_global_attr(text: str, name: str, value: str) -> str:
    pattern = re.compile(rf"(:{re.escape(name)}\s*=\s*)\".*?\"(\s*;)", re.DOTALL)
    if not pattern.search(text):
        raise ValueError(f"Global attribute {name!r} was not found")
    return pattern.sub(rf'\1"{value}"\2', text, count=1)


def pad_or_trim(value: str, width: int) -> str:
    if len(value) >= width:
        return value[:width]
    return value.ljust(width)


@dataclass
class ElasticCdf:
    path: Path
    text: str
    xs_data_tab: np.ndarray
    xs_data_base: list[int]
    xs_data_inc: list[int]
    xs_tab_index: list[list[int]]
    xs_min: list[list[float]]
    xs_max: list[list[float]]
    xs_spacing: list[list[str]]

    @classmethod
    def from_file(cls, path: str | Path) -> "ElasticCdf":
        cdf_path = Path(path)
        text = cdf_path.read_text(encoding="utf-8")
        xs_data_tab = np.array(parse_float_array(text, "xs_data_tab"), dtype=float)
        xs_data_base = parse_int_array(text, "xs_data_base")
        xs_data_inc = parse_int_array(text, "xs_data_inc")
        xs_tab_flat = parse_int_array(text, "xs_tab_index")
        xs_min_flat = parse_float_array(text, "xs_min")
        xs_max_flat = parse_float_array(text, "xs_max")
        xs_spacing_flat = parse_array(text, "xs_spacing")
        xs_tab_index = [
            xs_tab_flat[index * RANK_IND : (index + 1) * RANK_IND] for index in range(NUM_DEP_VAR)
        ]
        xs_min = [xs_min_flat[index * RANK_IND : (index + 1) * RANK_IND] for index in range(NUM_DEP_VAR)]
        xs_max = [xs_max_flat[index * RANK_IND : (index + 1) * RANK_IND] for index in range(NUM_DEP_VAR)]
        xs_spacing = [
            xs_spacing_flat[index * RANK_IND0 : (index + 1) * RANK_IND0]
            for index in range(NUM_DEP_VAR)
        ]
        return cls(
            path=cdf_path,
            text=text,
            xs_data_tab=xs_data_tab,
            xs_data_base=xs_data_base,
            xs_data_inc=xs_data_inc,
            xs_tab_index=xs_tab_index,
            xs_min=xs_min,
            xs_max=xs_max,
            xs_spacing=xs_spacing,
        )

    def get_block(self, index: int) -> np.ndarray:
        base = self.xs_data_base[index]
        size = self.xs_data_inc[index]
        return self.xs_data_tab[base : base + size].copy()

    def set_block(self, index: int, values: np.ndarray) -> None:
        base = self.xs_data_base[index]
        size = self.xs_data_inc[index]
        flat_values = np.asarray(values, dtype=float).ravel()
        if flat_values.size != size:
            raise ValueError(
                f"Block {index} expects {size} values, but received {flat_values.size}"
            )
        self.xs_data_tab[base : base + size] = flat_values

    def get_axis(self, block_index: int, dim_index: int) -> np.ndarray:
        count = self.xs_tab_index[block_index][dim_index]
        value_min = self.xs_min[block_index][dim_index]
        value_max = self.xs_max[block_index][dim_index]
        spacing = self.xs_spacing[block_index][dim_index]
        return make_axis(count, value_min, value_max, spacing)

    def render(
        self,
        *,
        xs_name: str | None = None,
        description: str | None = None,
        data_version: str | None = None,
    ) -> str:
        rendered = replace_array_block(self.text, "xs_data_tab", format_float_block(self.xs_data_tab))
        if xs_name is not None:
            rendered = replace_char_scalar(rendered, "xs_name", pad_or_trim(xs_name, 24))
        if description is not None:
            rendered = replace_global_attr(rendered, "description", description)
        if data_version is not None:
            rendered = replace_global_attr(rendered, "data_version", data_version)
        return rendered

    def write(
        self,
        output_path: str | Path,
        *,
        xs_name: str | None = None,
        description: str | None = None,
        data_version: str | None = None,
    ) -> None:
        Path(output_path).write_text(
            self.render(xs_name=xs_name, description=description, data_version=data_version),
            encoding="utf-8",
        )


def cross_section_axis(cdf: ElasticCdf) -> np.ndarray:
    return cdf.get_axis(0, 0)


def transport_axes(cdf: ElasticCdf) -> tuple[np.ndarray, np.ndarray]:
    return cdf.get_axis(1, 0), cdf.get_axis(1, 1)


def angle_axes(cdf: ElasticCdf) -> tuple[np.ndarray, np.ndarray]:
    return cdf.get_axis(5, 0), cdf.get_axis(5, 1)


def load_runtime_elastic_tables(
    path: str | Path,
    *,
    sigma_scale: float = 1.0,
) -> dict[str, np.ndarray]:
    cdf = ElasticCdf.from_file(path)
    sigma = cdf.get_block(0) * sigma_scale
    prob_axis, angle_energy = angle_axes(cdf)
    angle_grid = reshape_scattering_angles(cdf.get_block(5), cdf).T
    return {
        "energy_grid_sigma": cross_section_axis(cdf),
        "sigma_elastic": sigma,
        "prob_grid": prob_axis,
        "energy_grid_angle": angle_energy,
        "angle_cdf": angle_grid,
    }


def reshape_transport(flat_values: np.ndarray, cdf: ElasticCdf) -> np.ndarray:
    n_energy, n_temp, _ = cdf.xs_tab_index[1]
    return np.asarray(flat_values, dtype=float).reshape(n_temp, n_energy)


def reshape_scattering_angles(flat_values: np.ndarray, cdf: ElasticCdf) -> np.ndarray:
    n_prob, n_energy, _ = cdf.xs_tab_index[5]
    return np.asarray(flat_values, dtype=float).reshape(n_energy, n_prob)


def load_cross_section_csv(path: str | Path, value_column: str = "cross_section_cm2") -> np.ndarray:
    rows = list(csv.DictReader(Path(path).open("r", encoding="utf-8", newline="")))
    if not rows:
        raise ValueError(f"No rows found in {path}")
    return np.array([float(row[value_column]) for row in rows], dtype=float)


def load_scattering_angle_csv(path: str | Path) -> np.ndarray:
    rows = list(csv.DictReader(Path(path).open("r", encoding="utf-8", newline="")))
    if not rows:
        raise ValueError(f"No rows found in {path}")
    return np.array([float(row["angle_rad"]) for row in rows], dtype=float)


def load_transport_csv(path: str | Path) -> dict[str, np.ndarray]:
    rows = list(csv.DictReader(Path(path).open("r", encoding="utf-8", newline="")))
    if not rows:
        raise ValueError(f"No rows found in {path}")
    return {
        "reaction_rate": np.array([float(row["reaction_rate_cm3_s"]) for row in rows], dtype=float),
        "I_1_0": np.array([float(row["I_1_0_cm3_s"]) for row in rows], dtype=float),
        "I_1_1_up": np.array([float(row["I_1_1_up_cm4_s2"]) for row in rows], dtype=float),
        "I_1_2_up2": np.array([float(row["I_1_2_up2_cm5_s3"]) for row in rows], dtype=float),
    }


def load_scalars_json(path: str | Path) -> dict[str, float]:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def probability_axis(cdf: ElasticCdf) -> np.ndarray:
    return cdf.get_axis(5, 0)


def angle_momentum_transfer_factor(theta_values: np.ndarray, prob_axis: np.ndarray) -> float:
    return float(integrate.simpson(1.0 - np.cos(theta_values), x=prob_axis))


def reshape_runtime_scattering_angles(flat_values: np.ndarray, cdf: ElasticCdf) -> np.ndarray:
    return reshape_scattering_angles(flat_values, cdf)


def flatten_runtime_scattering_angles(angle_grid: np.ndarray) -> np.ndarray:
    return np.asarray(angle_grid, dtype=float).ravel()


def evaluate_angle_transport_ratio(scattering_angle: np.ndarray, cdf: ElasticCdf) -> np.ndarray:
    angle_grid = reshape_runtime_scattering_angles(scattering_angle, cdf)
    prob_axis = probability_axis(cdf)
    return np.array(
        [angle_momentum_transfer_factor(theta_values, prob_axis) for theta_values in angle_grid],
        dtype=float,
    )


def build_scaled_angle_table(
    cdf: ElasticCdf,
    reference_scattering_angle: np.ndarray,
    target_transport_ratio: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    ref_grid = reshape_runtime_scattering_angles(reference_scattering_angle, cdf)
    target = np.asarray(target_transport_ratio, dtype=float)
    if target.size != ref_grid.shape[0]:
        raise ValueError(
            f"Expected {ref_grid.shape[0]} target transport-ratio values, received {target.size}"
        )

    scaled_grid = np.zeros_like(ref_grid)
    scale_factors = np.zeros(ref_grid.shape[0], dtype=float)
    prob_axis = probability_axis(cdf)

    for index in range(ref_grid.shape[0]):
        theta_values = ref_grid[index]
        target_value = float(target[index])
        lo = 0.0
        hi = 4.0
        while angle_momentum_transfer_factor(np.clip(hi * theta_values, 0.0, np.pi), prob_axis) < target_value:
            hi *= 2.0
            if hi > 1024.0:
                raise ValueError(
                    f"Failed to bracket transport ratio {target_value:.6e} for angle slice {index}"
                )

        for _ in range(80):
            mid = 0.5 * (lo + hi)
            trial = np.clip(mid * theta_values, 0.0, np.pi)
            if angle_momentum_transfer_factor(trial, prob_axis) < target_value:
                lo = mid
            else:
                hi = mid

        scale_factors[index] = 0.5 * (lo + hi)
        scaled_grid[index] = np.clip(scale_factors[index] * theta_values, 0.0, np.pi)

    return flatten_runtime_scattering_angles(scaled_grid), scale_factors


def compute_transport_tables(
    cdf: ElasticCdf,
    cross_section: np.ndarray,
    scattering_angle: np.ndarray,
    *,
    xi_points: int = 1000,
) -> dict[str, np.ndarray | float]:
    sigma_energy = cross_section_axis(cdf)
    random_axis, angle_energy = angle_axes(cdf)
    angle_grid = reshape_scattering_angles(scattering_angle, cdf)
    r_theta = np.zeros(angle_grid.shape[0], dtype=float)
    for index, theta_values in enumerate(angle_grid):
        r_theta[index] = float(integrate.simpson(1.0 - np.cos(theta_values), x=random_axis))

    r_theta_interp = interp1d(
        angle_energy,
        r_theta,
        kind="cubic",
        bounds_error=False,
        fill_value=(float(r_theta[0]), float(r_theta[-1])),
    )

    return compute_transport_tables_from_sigma_momentum(
        cdf,
        cross_section,
        sigma_momentum=None,
        sigma_energy=sigma_energy,
        sigma_momentum_interp=r_theta_interp,
        xi_points=xi_points,
    )


def compute_transport_tables_from_sigma_momentum(
    cdf: ElasticCdf,
    cross_section: np.ndarray,
    sigma_momentum: np.ndarray | None,
    *,
    sigma_energy: np.ndarray | None = None,
    sigma_momentum_interp=None,
    xi_points: int = 1000,
) -> dict[str, np.ndarray | float]:
    if sigma_energy is None:
        sigma_energy = cross_section_axis(cdf)
    sigma_energy = np.asarray(sigma_energy, dtype=float)

    sigma_total_interp = interp1d(
        sigma_energy,
        np.asarray(cross_section, dtype=float),
        kind="cubic",
        bounds_error=False,
        fill_value=(float(cross_section[0]), float(cross_section[-1])),
    )
    if sigma_momentum_interp is None:
        if sigma_momentum is None:
            raise ValueError("sigma_momentum or sigma_momentum_interp must be provided")
        sigma_momentum = np.asarray(sigma_momentum, dtype=float)
        sigma_momentum_interp = interp1d(
            sigma_energy,
            sigma_momentum,
            kind="cubic",
            bounds_error=False,
            fill_value=(float(sigma_momentum[0]), float(sigma_momentum[-1])),
        )

    energy_axis, temp_axis = transport_axes(cdf)
    rr_table = np.zeros((temp_axis.size, energy_axis.size), dtype=float)
    i10_table = np.zeros_like(rr_table)
    i11_table = np.zeros_like(rr_table)
    i12_table = np.zeros_like(rr_table)

    for temp_index, temp_value in enumerate(temp_axis):
        a_beta = C0_CM_PER_S * np.sqrt(temp_value)
        for energy_index, energy_value in enumerate(energy_axis):
            v_alpha = C0_CM_PER_S * np.sqrt(energy_value)
            delta = np.sqrt(energy_value / temp_value)
            xi_max = delta + 6.0
            xi_grid = np.linspace(0.0, xi_max, xi_points)
            rel_energy = temp_value * xi_grid**2

            sigma_total = sigma_total_interp(rel_energy)
            sigma_mt = sigma_momentum_interp(rel_energy)

            exp_minus = np.exp(-(xi_grid - delta) ** 2)
            exp_plus = np.exp(-(xi_grid + delta) ** 2)
            diff_exp = exp_minus - exp_plus
            sum_exp = exp_minus + exp_plus

            prefactor = 1.0 / (np.sqrt(np.pi) * v_alpha)
            rr_table[temp_index, energy_index] = prefactor * (a_beta**2) * float(
                integrate.simpson(xi_grid**2 * sigma_total * diff_exp, x=xi_grid)
            )
            i10_table[temp_index, energy_index] = prefactor * (a_beta**2) * float(
                integrate.simpson(xi_grid**2 * sigma_mt * diff_exp, x=xi_grid)
            )
            i11_table[temp_index, energy_index] = prefactor * (a_beta**3) * float(
                integrate.simpson(xi_grid**3 * sigma_mt * sum_exp, x=xi_grid)
            )
            i12_table[temp_index, energy_index] = prefactor * (a_beta**4) * float(
                integrate.simpson(xi_grid**4 * sigma_mt * diff_exp, x=xi_grid)
            )

    return {
        "reaction_rate": rr_table.ravel(),
        "I_1_0": i10_table.ravel(),
        "I_1_1_up": i11_table.ravel(),
        "I_1_2_up2": i12_table.ravel(),
        "sigv_max": float(np.max(rr_table)),
    }


def apply_blocks(
    cdf: ElasticCdf,
    *,
    cross_section: np.ndarray,
    scattering_angle: np.ndarray,
    reaction_rate: np.ndarray | None = None,
    i_1_0: np.ndarray | None = None,
    i_1_1_up: np.ndarray | None = None,
    i_1_2_up2: np.ndarray | None = None,
    sigv_max: float | None = None,
    angle_min: float | None = None,
) -> None:
    cdf.set_block(0, cross_section)
    cdf.set_block(5, scattering_angle)
    if reaction_rate is not None:
        cdf.set_block(1, reaction_rate)
    if i_1_0 is not None:
        cdf.set_block(2, i_1_0)
    if i_1_1_up is not None:
        cdf.set_block(3, i_1_1_up)
    if i_1_2_up2 is not None:
        cdf.set_block(4, i_1_2_up2)
    if sigv_max is not None:
        cdf.set_block(6, np.array([sigv_max], dtype=float))
    if angle_min is not None:
        cdf.set_block(7, np.array([angle_min], dtype=float))
