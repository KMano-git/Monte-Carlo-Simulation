#!/usr/bin/env python3
"""Krstic D+ + D differential-cross-section utilities."""

from __future__ import annotations

import ast
import json
import math
import re
import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy import integrate


BOHR_AREA_CM2 = 2.80028e-17
DEFAULT_THETA_MIN_RAD = 1.0e-8
DEFAULT_THETA_MAX_RAD = math.pi - 1.0e-8
DEFAULT_THETA_SPLIT_RAD = 1.0e-2
DEFAULT_LOG_POINTS = 8192
DEFAULT_LINEAR_POINTS = 8192
DEFAULT_POLE_MARGIN = 1.0e-3
DEFAULT_NEGATIVE_REL_TOL = 1.0e-10
MANUAL_SECTION_MARKER = "## 11. DCS係数表の手動入力欄"


class NegativePureElasticError(RuntimeError):
    """Raised when the reconstructed pure-elastic kernel is significantly negative."""


@dataclass(frozen=True)
class ChannelCoefficients:
    a: tuple[float, ...]
    b: tuple[float, ...]
    A: float
    B: float = 0.0
    C: float = 0.0


@dataclass(frozen=True)
class EnergySlice:
    energy_cm_ev: float
    elastic: ChannelCoefficients
    spin_exchange: ChannelCoefficients


def _to_channel(entry: dict[str, object]) -> ChannelCoefficients:
    return ChannelCoefficients(
        a=tuple(float(value) for value in entry["a"]),
        b=tuple(float(value) for value in entry["b"]),
        A=float(entry["A"]),
        B=float(entry.get("B", 0.0)),
        C=float(entry.get("C", 0.0)),
    )


def parse_krstic_dcs_markdown(path: str | Path) -> dict[str, dict[str, dict[str, object]]]:
    markdown_path = Path(path)
    text = markdown_path.read_text(encoding="utf-8")
    if MANUAL_SECTION_MARKER not in text:
        raise ValueError(f"{MANUAL_SECTION_MARKER!r} was not found in {markdown_path}")

    body = text.split(MANUAL_SECTION_MARKER, maxsplit=1)[1]
    parsed: dict[str, dict[str, dict[str, object]]] = {}
    current_energy: str | None = None
    current_channel: str | None = None

    for line in body.splitlines():
        energy_match = re.match(r"^## E =\s*([0-9.]+)\s*eV", line)
        if energy_match:
            current_energy = f"{float(energy_match.group(1)):.4f}"
            parsed[current_energy] = {}
            current_channel = None
            continue

        channel_match = re.match(r"^###\s+(Elastic|Spin Exchange)\s*$", line)
        if channel_match:
            if current_energy is None:
                raise ValueError("Encountered channel header before energy header.")
            current_channel = (
                "elastic" if channel_match.group(1) == "Elastic" else "spin_exchange"
            )
            parsed[current_energy][current_channel] = {"a": [], "b": []}
            continue

        value_match = re.match(r"^\s*-\s+([A-Za-z_]+):\s*(.*)$", line)
        if value_match and current_energy is not None and current_channel is not None:
            key = value_match.group(1)
            raw_value = value_match.group(2).strip()
            channel_entry = parsed[current_energy][current_channel]
            if key.endswith("_values"):
                coeff_key = key[0]
                channel_entry[coeff_key] = [
                    float(value) for value in ast.literal_eval(raw_value.replace("E", "e"))
                ]
            elif key in {"A", "B", "C"}:
                channel_entry[key] = float(raw_value.replace("E", "e"))

    for energy_key, channels in parsed.items():
        for channel_name in ("elastic", "spin_exchange"):
            if channel_name not in channels:
                raise ValueError(f"Missing {channel_name} coefficients for E={energy_key} eV")
            required = {"a", "b", "A"}
            missing = required.difference(channels[channel_name])
            if missing:
                raise ValueError(
                    f"Missing keys {sorted(missing)} for {channel_name} at E={energy_key} eV"
                )
            channels[channel_name].setdefault("B", 0.0)
            channels[channel_name].setdefault("C", 0.0)
    return parsed


def write_coefficients_json(
    markdown_path: str | Path,
    output_path: str | Path,
) -> dict[str, dict[str, dict[str, object]]]:
    parsed = parse_krstic_dcs_markdown(markdown_path)
    Path(output_path).write_text(
        json.dumps(parsed, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return parsed


def load_coefficients_json(path: str | Path) -> list[EnergySlice]:
    raw = json.loads(Path(path).read_text(encoding="utf-8"))
    slices = [
        EnergySlice(
            energy_cm_ev=float(energy_key),
            elastic=_to_channel(energy_entry["elastic"]),
            spin_exchange=_to_channel(energy_entry["spin_exchange"]),
        )
        for energy_key, energy_entry in raw.items()
    ]
    slices.sort(key=lambda item: item.energy_cm_ev)
    return slices


def build_theta_grid(
    *,
    theta_min: float = DEFAULT_THETA_MIN_RAD,
    theta_split: float = DEFAULT_THETA_SPLIT_RAD,
    theta_max: float = DEFAULT_THETA_MAX_RAD,
    log_points: int = DEFAULT_LOG_POINTS,
    linear_points: int = DEFAULT_LINEAR_POINTS,
) -> np.ndarray:
    if not (0.0 < theta_min < theta_split < theta_max < math.pi):
        raise ValueError("Theta grid bounds must satisfy 0 < theta_min < theta_split < theta_max < pi.")
    log_grid = np.geomspace(theta_min, theta_split, log_points, endpoint=False)
    linear_grid = np.linspace(theta_split, theta_max, linear_points)
    return np.unique(np.concatenate((log_grid, linear_grid)))


def _series_in_log_theta(log_theta: np.ndarray, coeffs: tuple[float, ...]) -> np.ndarray:
    total = np.zeros_like(log_theta)
    for power, coeff in enumerate(coeffs):
        total += coeff * log_theta**power
    return total


def channel_pole_thetas(
    coeffs: ChannelCoefficients,
    *,
    theta_min: float = DEFAULT_THETA_MIN_RAD,
    theta_max: float = DEFAULT_THETA_MAX_RAD,
) -> np.ndarray:
    if not coeffs.b:
        return np.array([], dtype=float)
    polynomial = [1.0, *coeffs.b]
    roots = np.roots(list(reversed(polynomial)))
    real_roots = []
    for root in roots:
        if abs(root.imag) > 1.0e-10:
            continue
        theta_value = math.exp(float(root.real))
        if theta_min < theta_value < theta_max:
            real_roots.append(theta_value)
    return np.array(sorted(real_roots), dtype=float)


def evaluate_channel_kernel(
    theta_rad: np.ndarray,
    coeffs: ChannelCoefficients,
    *,
    exponent_clip: float = 700.0,
) -> np.ndarray:
    theta_values = np.asarray(theta_rad, dtype=float)
    if np.any(theta_values <= 0.0) or np.any(theta_values >= math.pi):
        raise ValueError("theta must satisfy 0 < theta < pi")

    log_theta = np.log(theta_values)
    numerator = _series_in_log_theta(log_theta, coeffs.a)
    denominator = np.ones_like(log_theta)
    for power, coeff in enumerate(coeffs.b, start=1):
        denominator += coeff * log_theta**power
    if np.any(np.abs(denominator) < 1.0e-12):
        raise ZeroDivisionError("Encountered near-singular Krstic DCS denominator.")

    angular_prefactor = (
        coeffs.A
        + coeffs.B * (1.0 - np.cos(theta_values))
        + coeffs.C * np.sin(theta_values) ** 2
    )
    exponent = np.clip(numerator / denominator, -exponent_clip, exponent_clip)
    kernel = angular_prefactor * np.exp(exponent)
    if not np.all(np.isfinite(kernel)):
        raise ValueError("Non-finite Krstic DCS kernel value encountered.")
    return kernel


def evaluate_differential_cross_section(
    theta_rad: np.ndarray,
    coeffs: ChannelCoefficients,
    *,
    exponent_clip: float = 700.0,
) -> np.ndarray:
    theta_values = np.asarray(theta_rad, dtype=float)
    return evaluate_channel_kernel(theta_values, coeffs, exponent_clip=exponent_clip) / (
        2.0 * np.pi * np.sin(theta_values)
    )


def _effective_pure_kernel(
    elastic_total: np.ndarray,
    spin_exchange: np.ndarray,
    *,
    energy_cm_ev: float,
    negative_policy: str,
    negative_rel_tol: float,
) -> tuple[np.ndarray, dict[str, float]]:
    raw_pure = elastic_total - spin_exchange
    scale = max(1.0, float(np.max(np.abs(elastic_total))), float(np.max(np.abs(spin_exchange))))
    allowed_negative = negative_rel_tol * scale
    negative_mask = raw_pure < 0.0
    min_negative = float(np.min(raw_pure))

    diagnostics = {
        "negative_min_raw_au": min_negative,
        "negative_rel_tol": float(negative_rel_tol),
        "clipped_negative_points": float(np.count_nonzero(negative_mask)),
    }
    if min_negative >= 0.0:
        return raw_pure, diagnostics

    if min_negative < -allowed_negative:
        message = (
            f"g_pure became significantly negative at E_cm={energy_cm_ev:.6g} eV: "
            f"min={min_negative:.6e} a0^2."
        )
        if negative_policy == "assert":
            raise NegativePureElasticError(message)
        if negative_policy == "warn-clip":
            warnings.warn(message, stacklevel=2)
        elif negative_policy != "clip":
            raise ValueError(
                "negative_policy must be one of 'assert', 'warn-clip', or 'clip'."
            )
    elif negative_policy == "warn-clip":
        warnings.warn(
            f"Clipping small negative g_pure values at E_cm={energy_cm_ev:.6g} eV "
            f"(min={min_negative:.6e} a0^2).",
            stacklevel=2,
        )

    clipped = np.maximum(raw_pure, 0.0)
    diagnostics["clipped_negative_area_proxy_au"] = float(np.sum(np.maximum(-raw_pure, 0.0)))
    return clipped, diagnostics


def compute_slice_observables(
    energy_slice: EnergySlice,
    prob_axis: np.ndarray,
    *,
    theta_grid: np.ndarray | None = None,
    theta_min: float = DEFAULT_THETA_MIN_RAD,
    theta_max: float = DEFAULT_THETA_MAX_RAD,
    theta_split: float = DEFAULT_THETA_SPLIT_RAD,
    log_points: int = DEFAULT_LOG_POINTS,
    linear_points: int = DEFAULT_LINEAR_POINTS,
    pole_margin: float = DEFAULT_POLE_MARGIN,
    negative_policy: str = "assert",
    negative_rel_tol: float = DEFAULT_NEGATIVE_REL_TOL,
    exponent_clip: float = 700.0,
) -> dict[str, object]:
    pole_thetas = np.concatenate(
        (
            channel_pole_thetas(energy_slice.elastic, theta_min=theta_min, theta_max=theta_max),
            channel_pole_thetas(
                energy_slice.spin_exchange,
                theta_min=theta_min,
                theta_max=theta_max,
            ),
        )
    )
    theta_min_used = float(theta_min)
    if pole_thetas.size:
        theta_min_used = max(theta_min_used, float(np.max(pole_thetas) * (1.0 + pole_margin)))
    if theta_min_used >= theta_max:
        raise ValueError(
            f"No valid theta interval remains at E_cm={energy_slice.energy_cm_ev:.6g} eV."
        )

    if theta_grid is None:
        theta_values = build_theta_grid(
            theta_min=theta_min_used,
            theta_split=max(theta_split, theta_min_used * 1.01),
            theta_max=theta_max,
            log_points=log_points,
            linear_points=linear_points,
        )
    else:
        theta_values = np.asarray(theta_grid, dtype=float)
        theta_values = theta_values[(theta_values >= theta_min_used) & (theta_values <= theta_max)]
        if theta_values.size < 4:
            raise ValueError(
                f"Theta grid became too short after pole filtering at E_cm={energy_slice.energy_cm_ev:.6g} eV."
            )

    g_elastic_total = evaluate_channel_kernel(
        theta_values,
        energy_slice.elastic,
        exponent_clip=exponent_clip,
    )
    g_spin_exchange = evaluate_channel_kernel(
        theta_values,
        energy_slice.spin_exchange,
        exponent_clip=exponent_clip,
    )
    g_pure, negative_diag = _effective_pure_kernel(
        g_elastic_total,
        g_spin_exchange,
        energy_cm_ev=energy_slice.energy_cm_ev,
        negative_policy=negative_policy,
        negative_rel_tol=negative_rel_tol,
    )

    sigma_t_elastic_total_au = float(integrate.simpson(g_elastic_total, x=theta_values))
    sigma_t_spin_exchange_au = float(integrate.simpson(g_spin_exchange, x=theta_values))
    sigma_mt_elastic_total_au = float(
        integrate.simpson((1.0 - np.cos(theta_values)) * g_elastic_total, x=theta_values)
    )
    sigma_vi_elastic_total_au = float(
        integrate.simpson(np.sin(theta_values) ** 2 * g_elastic_total, x=theta_values)
    )
    sigma_t_pure_au = float(integrate.simpson(g_pure, x=theta_values))
    sigma_mt_pure_au = float(
        integrate.simpson((1.0 - np.cos(theta_values)) * g_pure, x=theta_values)
    )
    sigma_vi_pure_au = float(
        integrate.simpson(np.sin(theta_values) ** 2 * g_pure, x=theta_values)
    )
    if sigma_t_pure_au <= 0.0:
        raise ValueError(
            f"Non-positive pure-elastic total cross section at E_cm={energy_slice.energy_cm_ev:.6g} eV."
        )

    cumulative = integrate.cumulative_trapezoid(g_pure, theta_values, initial=0.0)
    cumulative = np.maximum.accumulate(cumulative)
    cumulative /= cumulative[-1]

    unique_cdf, unique_indices = np.unique(cumulative, return_index=True)
    theta_support = theta_values[unique_indices]
    if unique_cdf[0] > 0.0:
        unique_cdf = np.insert(unique_cdf, 0, 0.0)
        theta_support = np.insert(theta_support, 0, theta_min_used)
    if unique_cdf[-1] < 1.0:
        unique_cdf = np.append(unique_cdf, 1.0)
        theta_support = np.append(theta_support, theta_values[-1])

    inverse_cdf = np.interp(prob_axis, unique_cdf, theta_support)
    inverse_cdf[0] = 0.0
    inverse_cdf[-1] = math.pi
    transport_ratio_from_cdf = float(
        integrate.simpson(1.0 - np.cos(inverse_cdf), x=np.asarray(prob_axis, dtype=float))
    )

    return {
        "energy_cm_ev": energy_slice.energy_cm_ev,
        "energy_lab_ev": 2.0 * energy_slice.energy_cm_ev,
        "theta_min_used_rad": theta_min_used,
        "pole_thetas_rad": pole_thetas.tolist(),
        "sigma_t_elastic_total_au": sigma_t_elastic_total_au,
        "sigma_t_spin_exchange_au": sigma_t_spin_exchange_au,
        "sigma_mt_elastic_total_au": sigma_mt_elastic_total_au,
        "sigma_vi_elastic_total_au": sigma_vi_elastic_total_au,
        "sigma_t_pure_au": sigma_t_pure_au,
        "sigma_mt_pure_au": sigma_mt_pure_au,
        "sigma_vi_pure_au": sigma_vi_pure_au,
        "sigma_t_elastic_total_cm2": sigma_t_elastic_total_au * BOHR_AREA_CM2,
        "sigma_t_spin_exchange_cm2": sigma_t_spin_exchange_au * BOHR_AREA_CM2,
        "sigma_mt_elastic_total_cm2": sigma_mt_elastic_total_au * BOHR_AREA_CM2,
        "sigma_vi_elastic_total_cm2": sigma_vi_elastic_total_au * BOHR_AREA_CM2,
        "sigma_t_pure_cm2": sigma_t_pure_au * BOHR_AREA_CM2,
        "sigma_mt_pure_cm2": sigma_mt_pure_au * BOHR_AREA_CM2,
        "sigma_vi_pure_cm2": sigma_vi_pure_au * BOHR_AREA_CM2,
        "transport_ratio_pure": sigma_mt_pure_au / sigma_t_pure_au,
        "transport_ratio_from_cdf": transport_ratio_from_cdf,
        "inverse_cdf_rad": inverse_cdf,
        "theta_grid_rad": theta_values,
        "g_elastic_total_au": g_elastic_total,
        "g_spin_exchange_au": g_spin_exchange,
        "g_pure_au": g_pure,
        **negative_diag,
    }


def build_tabulated_observables(
    energy_slices: list[EnergySlice],
    prob_axis: np.ndarray,
    *,
    theta_grid: np.ndarray | None = None,
    theta_min: float = DEFAULT_THETA_MIN_RAD,
    theta_max: float = DEFAULT_THETA_MAX_RAD,
    theta_split: float = DEFAULT_THETA_SPLIT_RAD,
    log_points: int = DEFAULT_LOG_POINTS,
    linear_points: int = DEFAULT_LINEAR_POINTS,
    pole_margin: float = DEFAULT_POLE_MARGIN,
    negative_policy: str = "assert",
    negative_rel_tol: float = DEFAULT_NEGATIVE_REL_TOL,
    exponent_clip: float = 700.0,
) -> dict[str, object]:
    rows = [
        compute_slice_observables(
            energy_slice,
            prob_axis,
            theta_grid=theta_grid,
            theta_min=theta_min,
            theta_max=theta_max,
            theta_split=theta_split,
            log_points=log_points,
            linear_points=linear_points,
            pole_margin=pole_margin,
            negative_policy=negative_policy,
            negative_rel_tol=negative_rel_tol,
            exponent_clip=exponent_clip,
        )
        for energy_slice in energy_slices
    ]
    return {
        "rows": rows,
        "energy_cm_ev": np.array([row["energy_cm_ev"] for row in rows], dtype=float),
        "energy_lab_ev": np.array([row["energy_lab_ev"] for row in rows], dtype=float),
        "sigma_t_pure_cm2": np.array([row["sigma_t_pure_cm2"] for row in rows], dtype=float),
        "sigma_mt_pure_cm2": np.array([row["sigma_mt_pure_cm2"] for row in rows], dtype=float),
        "sigma_vi_pure_cm2": np.array([row["sigma_vi_pure_cm2"] for row in rows], dtype=float),
        "sigma_t_elastic_total_cm2": np.array(
            [row["sigma_t_elastic_total_cm2"] for row in rows],
            dtype=float,
        ),
        "sigma_mt_elastic_total_cm2": np.array(
            [row["sigma_mt_elastic_total_cm2"] for row in rows],
            dtype=float,
        ),
        "sigma_vi_elastic_total_cm2": np.array(
            [row["sigma_vi_elastic_total_cm2"] for row in rows],
            dtype=float,
        ),
        "sigma_t_spin_exchange_cm2": np.array(
            [row["sigma_t_spin_exchange_cm2"] for row in rows],
            dtype=float,
        ),
        "inverse_cdf_rad": np.vstack([row["inverse_cdf_rad"] for row in rows]),
    }


def interpolate_inverse_cdf(
    tabulated_energy_cm_ev: np.ndarray,
    tabulated_inverse_cdf: np.ndarray,
    target_energy_cm_ev: np.ndarray,
) -> np.ndarray:
    energy_in = np.asarray(tabulated_energy_cm_ev, dtype=float)
    inverse_in = np.asarray(tabulated_inverse_cdf, dtype=float)
    energy_out = np.asarray(target_energy_cm_ev, dtype=float)
    if inverse_in.shape[0] != energy_in.size:
        raise ValueError("Inverse-CDF rows must match the tabulated energy axis.")

    log_energy_in = np.log(energy_in)
    log_energy_out = np.log(np.clip(energy_out, energy_in[0], energy_in[-1]))
    inverse_out = np.zeros((energy_out.size, inverse_in.shape[1]), dtype=float)
    for prob_index in range(inverse_in.shape[1]):
        inverse_out[:, prob_index] = np.interp(
            log_energy_out,
            log_energy_in,
            inverse_in[:, prob_index],
            left=float(inverse_in[0, prob_index]),
            right=float(inverse_in[-1, prob_index]),
        )
    inverse_out[:, 0] = 0.0
    inverse_out[:, -1] = math.pi
    return inverse_out


def interpolate_positive_observable(
    energy_in: np.ndarray,
    values_in: np.ndarray,
    energy_out: np.ndarray,
) -> np.ndarray:
    energy_src = np.asarray(energy_in, dtype=float)
    values_src = np.asarray(values_in, dtype=float)
    energy_dst = np.asarray(energy_out, dtype=float)
    if np.any(energy_src <= 0.0) or np.any(values_src <= 0.0):
        raise ValueError("Positive log-log interpolation requires strictly positive inputs.")

    log_energy_src = np.log(energy_src)
    log_values_src = np.log(values_src)
    log_energy_dst = np.log(np.clip(energy_dst, energy_src[0], energy_src[-1]))
    return np.exp(
        np.interp(
            log_energy_dst,
            log_energy_src,
            log_values_src,
            left=float(log_values_src[0]),
            right=float(log_values_src[-1]),
        )
    )


def interpolate_observables_to_lab_grid(
    tabulated: dict[str, object],
    energy_lab_ev: np.ndarray,
) -> dict[str, np.ndarray]:
    energy_cm_out = 0.5 * np.asarray(energy_lab_ev, dtype=float)
    energy_cm_in = np.asarray(tabulated["energy_cm_ev"], dtype=float)
    return {
        "sigma_t_pure_cm2": interpolate_positive_observable(
            energy_cm_in,
            np.asarray(tabulated["sigma_t_pure_cm2"], dtype=float),
            energy_cm_out,
        ),
        "sigma_mt_pure_cm2": interpolate_positive_observable(
            energy_cm_in,
            np.asarray(tabulated["sigma_mt_pure_cm2"], dtype=float),
            energy_cm_out,
        ),
        "sigma_vi_pure_cm2": interpolate_positive_observable(
            energy_cm_in,
            np.asarray(tabulated["sigma_vi_pure_cm2"], dtype=float),
            energy_cm_out,
        ),
        "inverse_cdf_rad": interpolate_inverse_cdf(
            energy_cm_in,
            np.asarray(tabulated["inverse_cdf_rad"], dtype=float),
            energy_cm_out,
        ),
    }


def load_or_build_coefficients(
    *,
    coeff_json_path: str | Path,
    markdown_path: str | Path | None = None,
) -> list[EnergySlice]:
    coeff_path = Path(coeff_json_path)
    if not coeff_path.exists():
        if markdown_path is None:
            raise FileNotFoundError(coeff_path)
        write_coefficients_json(markdown_path, coeff_path)
    return load_coefficients_json(coeff_path)


def main() -> None:
    this_dir = Path(__file__).resolve().parent
    markdown_path = this_dir / "DpD_fit_memo_v2.md"
    coeff_json_path = this_dir / "Krstic" / "krstic_pure_dcs_coeffs.json"
    coeff_json_path.parent.mkdir(parents=True, exist_ok=True)
    write_coefficients_json(markdown_path, coeff_json_path)
    print(f"Wrote {coeff_json_path}")


if __name__ == "__main__":
    main()
