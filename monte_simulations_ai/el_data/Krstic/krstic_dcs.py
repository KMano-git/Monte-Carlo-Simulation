#!/usr/bin/env python3
"""Parse and evaluate Krstic D + D+ differential cross-section fits."""

from __future__ import annotations

import ast
import math
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy import integrate


BOHR_AREA_CM2 = 2.80028e-17
MANUAL_SECTION_MARKER = "## 11. DCS係数表の手動入力欄"


@dataclass(frozen=True)
class DifferentialFit:
    energy_ev: float
    channel: str
    a_label: str
    a_values: tuple[float, ...]
    b_label: str
    b_values: tuple[float, ...]
    A: float
    B: float
    C: float

    def theta_kernel_au(self, theta: np.ndarray) -> np.ndarray:
        """Return 2*pi*sin(theta)*d sigma / d Omega in a0^2."""
        theta_values = np.asarray(theta, dtype=float)
        if np.any(theta_values <= 0.0):
            raise ValueError("Scattering angle must be strictly positive.")

        log_theta = np.log(theta_values)
        numerator = np.zeros_like(theta_values)
        for power, coeff in enumerate(self.a_values):
            numerator += coeff * log_theta**power

        denominator = np.ones_like(theta_values)
        for power, coeff in enumerate(self.b_values, start=1):
            denominator += coeff * log_theta**power
        if np.any(np.abs(denominator) < 1.0e-12):
            raise ZeroDivisionError(
                f"Near-singular differential-fit denominator for {self.channel} at "
                f"E={self.energy_ev:.6g} eV."
            )

        angular_prefactor = (
            self.A
            + self.B * (1.0 - np.cos(theta_values))
            + self.C * np.sin(theta_values) ** 2
        )
        kernel = angular_prefactor * np.exp(numerator / denominator)

        if not np.all(np.isfinite(kernel)):
            raise ValueError(
                f"Non-finite differential-fit values for {self.channel} at "
                f"E={self.energy_ev:.6g} eV."
            )
        if np.min(kernel) < 0.0:
            raise ValueError(
                f"Negative differential-fit values for {self.channel} at "
                f"E={self.energy_ev:.6g} eV."
            )
        return kernel

    def differential_cross_section_au_per_sr(self, theta: np.ndarray) -> np.ndarray:
        theta_values = np.asarray(theta, dtype=float)
        kernel = self.theta_kernel_au(theta_values)
        return kernel / (2.0 * np.pi * np.sin(theta_values))


def parse_krstic_dcs_markdown(path: str | Path) -> dict[str, list[DifferentialFit]]:
    markdown_path = Path(path)
    text = markdown_path.read_text(encoding="utf-8")
    if MANUAL_SECTION_MARKER not in text:
        raise ValueError(f"{MANUAL_SECTION_MARKER!r} was not found in {markdown_path}")

    body = text.split(MANUAL_SECTION_MARKER, maxsplit=1)[1]
    body = re.split(r"\n##\s+(?!E\s*=)", body, maxsplit=1)[0]
    fits_by_channel: dict[str, list[DifferentialFit]] = {
        "Elastic": [],
        "Spin Exchange": [],
    }

    current_energy: float | None = None
    current_channel: str | None = None
    current_values: dict[str, str] = {}

    def flush_current() -> None:
        nonlocal current_values
        if current_energy is None or current_channel is None:
            current_values = {}
            return

        required_keys = ("a_label", "a_values", "b_label", "b_values", "A", "B", "C")
        missing = [key for key in required_keys if key not in current_values]
        if missing:
            raise ValueError(
                f"Missing DCS keys for {current_channel} at E={current_energy:.6g} eV: {missing}"
            )

        fit = DifferentialFit(
            energy_ev=float(current_energy),
            channel=current_channel,
            a_label=current_values["a_label"],
            a_values=tuple(float(value) for value in ast.literal_eval(current_values["a_values"])),
            b_label=current_values["b_label"],
            b_values=tuple(float(value) for value in ast.literal_eval(current_values["b_values"])),
            A=float(current_values["A"]),
            B=float(current_values["B"]),
            C=float(current_values["C"]),
        )
        fits_by_channel[current_channel].append(fit)
        current_values = {}

    for line in body.splitlines():
        energy_match = re.match(r"^## E =\s*([0-9.]+)\s*eV", line)
        if energy_match:
            flush_current()
            current_energy = float(energy_match.group(1))
            current_channel = None
            continue

        channel_match = re.match(r"^###\s+(Elastic|Spin Exchange)\s*$", line)
        if channel_match:
            flush_current()
            current_channel = channel_match.group(1)
            continue

        value_match = re.match(r"^\s*-\s+([A-Za-z_]+):\s*(.*)$", line)
        if value_match and current_channel is not None:
            current_values[value_match.group(1)] = value_match.group(2).strip()

    flush_current()
    for channel, fits in fits_by_channel.items():
        fits.sort(key=lambda fit: fit.energy_ev)
        if not fits:
            raise ValueError(f"No DCS coefficients were found for channel {channel!r}")
    return fits_by_channel


def serialize_fits(
    fits_by_channel: dict[str, list[DifferentialFit]],
    *,
    source_markdown: str,
) -> dict[str, object]:
    return {
        "source_markdown": source_markdown,
        "formula": (
            "2*pi*sin(theta) d sigma / d Omega = "
            "[A + B(1-cos theta) + C sin^2 theta] * "
            "exp[(sum_i a_i (ln theta)^i) / (1 + sum_j b_j (ln theta)^j)]"
        ),
        "theta_units": "rad",
        "energy_units": "eV (center of mass / eV per amu for D + D+)",
        "cross_section_units": "a0^2",
        "channels": {
            channel.lower().replace(" ", "_"): [
                {
                    "energy_eV_amu": fit.energy_ev,
                    "a_label": fit.a_label,
                    "a_values": list(fit.a_values),
                    "b_label": fit.b_label,
                    "b_values": list(fit.b_values),
                    "A": fit.A,
                    "B": fit.B,
                    "C": fit.C,
                }
                for fit in fits
            ]
            for channel, fits in fits_by_channel.items()
        },
    }


def build_inverse_cdf(
    fit: DifferentialFit,
    prob_axis: np.ndarray,
    *,
    theta_floor: float = 1.0e-8,
    theta_points: int = 16384,
) -> dict[str, np.ndarray | float]:
    theta_values = np.geomspace(theta_floor, math.pi, theta_points)
    kernel = fit.theta_kernel_au(theta_values)
    return build_inverse_cdf_from_kernel(
        theta_values,
        kernel,
        prob_axis,
        energy_ev=fit.energy_ev,
        channel=fit.channel,
    )


def build_inverse_cdf_from_kernel(
    theta_values: np.ndarray,
    kernel_au: np.ndarray,
    prob_axis: np.ndarray,
    *,
    energy_ev: float,
    channel: str,
) -> dict[str, np.ndarray | float]:
    theta_values = np.asarray(theta_values, dtype=float)
    kernel = np.asarray(kernel_au, dtype=float)
    prob = np.asarray(prob_axis, dtype=float)
    if theta_values.ndim != 1 or kernel.ndim != 1:
        raise ValueError("Theta grid and kernel must be 1D arrays.")
    if theta_values.size != kernel.size:
        raise ValueError("Theta grid and kernel lengths differ.")
    if np.any(theta_values <= 0.0) or np.any(np.diff(theta_values) <= 0.0):
        raise ValueError("Theta grid must be strictly increasing and positive.")
    if np.any(kernel < 0.0) or not np.all(np.isfinite(kernel)):
        raise ValueError(
            f"Invalid differential-fit kernel for {channel} at E={energy_ev:.6g} eV."
        )

    cumulative = integrate.cumulative_trapezoid(kernel, theta_values, initial=0.0)
    sigma_total_au = float(cumulative[-1])
    if sigma_total_au <= 0.0:
        raise ValueError(
            f"Non-positive integrated cross section for {channel} at E={energy_ev:.6g} eV."
        )

    sigma_momentum_au = float(
        integrate.trapezoid((1.0 - np.cos(theta_values)) * kernel, x=theta_values)
    )
    probability_support = np.concatenate(([0.0], cumulative / sigma_total_au))
    theta_support = np.concatenate(([0.0], theta_values))
    inverse_cdf = np.interp(prob, probability_support, theta_support)
    inverse_cdf[0] = 0.0
    inverse_cdf[-1] = math.pi
    transport_ratio_from_inverse_cdf = float(
        integrate.simpson(1.0 - np.cos(inverse_cdf), x=prob)
    )

    return {
        "inverse_cdf_rad": inverse_cdf,
        "sigma_total_au": sigma_total_au,
        "sigma_total_cm2": sigma_total_au * BOHR_AREA_CM2,
        "sigma_momentum_au": sigma_momentum_au,
        "sigma_momentum_cm2": sigma_momentum_au * BOHR_AREA_CM2,
        "transport_ratio": sigma_momentum_au / sigma_total_au,
        "transport_ratio_from_inverse_cdf": transport_ratio_from_inverse_cdf,
    }


def build_tabulated_inverse_cdf(
    fits: list[DifferentialFit],
    prob_axis: np.ndarray,
    *,
    theta_floor: float = 1.0e-8,
    theta_points: int = 16384,
) -> tuple[np.ndarray, np.ndarray, list[dict[str, float]]]:
    energies = np.array([fit.energy_ev for fit in fits], dtype=float)
    inverse_cdf_table = np.zeros((len(fits), prob_axis.size), dtype=float)
    summaries: list[dict[str, float]] = []

    for index, fit in enumerate(fits):
        built = build_inverse_cdf(
            fit,
            prob_axis,
            theta_floor=theta_floor,
            theta_points=theta_points,
        )
        inverse_cdf_table[index] = np.asarray(built["inverse_cdf_rad"], dtype=float)
        summaries.append(
            {
                "energy_eV_amu": fit.energy_ev,
                "sigma_total_au": float(built["sigma_total_au"]),
                "sigma_total_cm2": float(built["sigma_total_cm2"]),
                "sigma_momentum_au": float(built["sigma_momentum_au"]),
                "sigma_momentum_cm2": float(built["sigma_momentum_cm2"]),
                "transport_ratio": float(built["transport_ratio"]),
            }
        )

    return energies, inverse_cdf_table, summaries


def interpolate_inverse_cdf_by_energy(
    tabulated_energy: np.ndarray,
    tabulated_inverse_cdf: np.ndarray,
    runtime_energy: np.ndarray,
) -> np.ndarray:
    energy_in = np.asarray(tabulated_energy, dtype=float)
    inverse_in = np.asarray(tabulated_inverse_cdf, dtype=float)
    energy_out = np.asarray(runtime_energy, dtype=float)

    if inverse_in.shape[0] != energy_in.size:
        raise ValueError("Energy axis length does not match inverse-CDF table rows.")

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


def interpolate_theta_kernel_by_energy(
    fits: list[DifferentialFit],
    theta_values: np.ndarray,
    runtime_energy: np.ndarray,
) -> np.ndarray:
    """Return DCS-first theta kernels interpolated in log(E)-log(kernel)."""
    if not fits:
        raise ValueError("At least one DCS fit is required.")
    sorted_fits = sorted(fits, key=lambda fit: fit.energy_ev)
    energy_in = np.array([fit.energy_ev for fit in sorted_fits], dtype=float)
    energy_out = np.asarray(runtime_energy, dtype=float)
    theta_grid = np.asarray(theta_values, dtype=float)
    if np.any(energy_in <= 0.0) or np.any(energy_out <= 0.0):
        raise ValueError("DCS energy interpolation requires positive energies.")

    kernel_in = np.vstack([fit.theta_kernel_au(theta_grid) for fit in sorted_fits])
    if np.any(kernel_in <= 0.0):
        raise ValueError("Log-kernel interpolation requires strictly positive DCS kernels.")

    log_energy_in = np.log(energy_in)
    log_kernel_in = np.log(kernel_in)
    log_energy_out = np.log(np.clip(energy_out, energy_in[0], energy_in[-1]))
    kernel_out = np.zeros((energy_out.size, theta_grid.size), dtype=float)

    for theta_index in range(theta_grid.size):
        kernel_out[:, theta_index] = np.exp(
            np.interp(
                log_energy_out,
                log_energy_in,
                log_kernel_in[:, theta_index],
                left=float(log_kernel_in[0, theta_index]),
                right=float(log_kernel_in[-1, theta_index]),
            )
        )
    return kernel_out


def build_dcs_first_inverse_cdf_by_energy(
    fits: list[DifferentialFit],
    prob_axis: np.ndarray,
    runtime_energy: np.ndarray,
    *,
    theta_floor: float = 1.0e-8,
    theta_points: int = 16384,
) -> tuple[np.ndarray, list[dict[str, float]]]:
    """Build inverse CDFs after interpolating p(theta, E) on the theta grid."""
    theta_values = np.geomspace(theta_floor, math.pi, theta_points)
    runtime_energy_values = np.asarray(runtime_energy, dtype=float)
    kernels = interpolate_theta_kernel_by_energy(fits, theta_values, runtime_energy_values)
    inverse_cdf_table = np.zeros(
        (runtime_energy_values.size, np.asarray(prob_axis).size),
        dtype=float,
    )
    summaries: list[dict[str, float]] = []

    for index, (energy_ev, kernel) in enumerate(
        zip(runtime_energy_values, kernels, strict=True)
    ):
        built = build_inverse_cdf_from_kernel(
            theta_values,
            kernel,
            prob_axis,
            energy_ev=float(energy_ev),
            channel="Elastic",
        )
        inverse_cdf_table[index] = np.asarray(built["inverse_cdf_rad"], dtype=float)
        summaries.append(
            {
                "energy_eV_amu": float(energy_ev),
                "sigma_total_au": float(built["sigma_total_au"]),
                "sigma_total_cm2": float(built["sigma_total_cm2"]),
                "sigma_momentum_au": float(built["sigma_momentum_au"]),
                "sigma_momentum_cm2": float(built["sigma_momentum_cm2"]),
                "transport_ratio": float(built["transport_ratio"]),
                "transport_ratio_from_inverse_cdf": float(
                    built["transport_ratio_from_inverse_cdf"]
                ),
            }
        )

    return inverse_cdf_table, summaries
