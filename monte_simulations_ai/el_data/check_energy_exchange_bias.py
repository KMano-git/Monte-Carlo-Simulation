#!/usr/bin/env python3
"""Check elastic energy-exchange bias from dd_00_elastic-compatible CDF files.

The primary diagnostic is the Maxwellian expectation of the table-lookup
elastic energy gain coefficient for neutral D:

    <sigma v DeltaE_n>  [eV cm^3/s]

For equal neutral/ion temperatures and zero drift this expectation should be
close to zero.  A nonzero value indicates bias from the I-kernel tables,
interpolation, or the runtime lookup formula.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy import integrate
from scipy.interpolate import RegularGridInterpolator, interp1d

from cdf_compat import (
    ElasticCdf,
    angle_axes,
    angle_momentum_transfer_factor,
    cross_section_axis,
    parse_float_array,
    probability_axis,
    reshape_scattering_angles,
    reshape_transport,
    transport_axes,
)


EV_TO_J = 1.602176634e-19
M_U_KG = 1.660538921e-27
M_D_KG = 2.01410177812 * M_U_KG
SPEED_CM_PER_M = 100.0

DEFAULT_CDFS = (
    "Bachmann/bachmann_dd_from_split_tables_compat.cdf",
    "dd_00_elastic_pure_el_angle_fixed.cdf",
    "Krstic/krstic_dd_total_elastic_integral_priority.cdf",
    "Krstic/krstic_dd_pure_dcs_compat.cdf",
)
DEFAULT_TEMPERATURES_EV = (0.5, 1.0, 1.8, 2.0, 3.0, 5.0, 10.0)


@dataclass
class TableInterpolator:
    energy_axis: np.ndarray
    temp_axis: np.ndarray
    interpolator: RegularGridInterpolator
    log_values: bool

    def __call__(self, energy: np.ndarray, temp: np.ndarray | float) -> np.ndarray:
        energy_array = np.asarray(energy, dtype=float)
        temp_array = np.asarray(temp, dtype=float)
        if temp_array.shape == ():
            temp_array = np.full_like(energy_array, float(temp_array))
        else:
            temp_array = np.broadcast_to(temp_array, energy_array.shape)

        energy_clipped = np.clip(energy_array, self.energy_axis[0], self.energy_axis[-1])
        temp_clipped = np.clip(temp_array, self.temp_axis[0], self.temp_axis[-1])
        points = np.column_stack(
            [np.log(energy_clipped).ravel(), np.log(temp_clipped).ravel()]
        )
        values = self.interpolator(points).reshape(energy_array.shape)
        if self.log_values:
            return np.exp(values)
        return values


def parse_temperature_list(values: list[str]) -> list[float]:
    temperatures: list[float] = []
    for value in values:
        for part in value.split(","):
            stripped = part.strip()
            if stripped:
                temperatures.append(float(stripped))
    if not temperatures:
        raise ValueError("No temperatures were provided.")
    return temperatures


def parse_pair_list(values: list[str]) -> list[tuple[float, float]]:
    pairs: list[tuple[float, float]] = []
    for value in values:
        for part in value.split():
            left, sep, right = part.partition(",")
            if not sep:
                raise ValueError(f"Temperature pair must be T_neutral,T_ion: {part!r}")
            pairs.append((float(left), float(right)))
    if not pairs:
        raise ValueError("No temperature pairs were provided.")
    return pairs


def block_multiplier(cdf: ElasticCdf, block_index: int) -> float:
    multipliers = parse_float_array(cdf.text, "xs_mult")
    offset = 4 * block_index
    if offset >= len(multipliers):
        raise ValueError(f"xs_mult does not contain block {block_index}")
    return float(multipliers[offset])


def build_table_interpolator(cdf: ElasticCdf, block_index: int) -> TableInterpolator:
    energy_axis, temp_axis = transport_axes(cdf)
    table_temp_energy = reshape_transport(cdf.get_block(block_index), cdf)
    table_energy_temp = table_temp_energy.T
    log_values = bool(np.all(table_energy_temp > 0.0))
    values = np.log(table_energy_temp) if log_values else table_energy_temp
    interpolator = RegularGridInterpolator(
        (np.log(energy_axis), np.log(temp_axis)),
        values,
        method="linear",
        bounds_error=False,
        fill_value=None,
    )
    return TableInterpolator(
        energy_axis=np.asarray(energy_axis, dtype=float),
        temp_axis=np.asarray(temp_axis, dtype=float),
        interpolator=interpolator,
        log_values=log_values,
    )


def maxwell_speed_pdf(speed_m_s: np.ndarray, temperature_ev: float) -> np.ndarray:
    sigma = np.sqrt(temperature_ev * EV_TO_J / M_D_KG)
    return (
        np.sqrt(2.0 / np.pi)
        * speed_m_s**2
        / sigma**3
        * np.exp(-(speed_m_s**2) / (2.0 * sigma**2))
    )


def speed_grid(temperature_ev: float, points: int, sigma_extent: float) -> np.ndarray:
    sigma = np.sqrt(temperature_ev * EV_TO_J / M_D_KG)
    return np.linspace(0.0, sigma_extent * sigma, points)


class EnergyExchangeDiagnostic:
    def __init__(self, cdf_path: Path) -> None:
        self.path = cdf_path
        self.cdf = ElasticCdf.from_file(cdf_path)
        self.reaction_rate = build_table_interpolator(self.cdf, 1)
        self.i_1_0 = build_table_interpolator(self.cdf, 2)
        self.i_1_1_up = build_table_interpolator(self.cdf, 3)
        self.i_1_2_up2 = build_table_interpolator(self.cdf, 4)
        self.reaction_rate_mult = block_multiplier(self.cdf, 1)
        self.i_1_0_mult = block_multiplier(self.cdf, 2)
        self.i_1_1_mult = block_multiplier(self.cdf, 3)
        self.i_1_2_mult = block_multiplier(self.cdf, 4)
        self.sigma_mt_from_angle = self._build_angle_sigma_mt_interpolator()

    def _build_angle_sigma_mt_interpolator(self):
        sigma_energy = cross_section_axis(self.cdf)
        sigma_total = self.cdf.get_block(0)
        prob_axis = probability_axis(self.cdf)
        _, angle_energy = angle_axes(self.cdf)
        angle_rows = reshape_scattering_angles(self.cdf.get_block(5), self.cdf)
        r_theta = np.array(
            [
                angle_momentum_transfer_factor(theta_values, prob_axis)
                for theta_values in angle_rows
            ],
            dtype=float,
        )

        sigma_interp = interp1d(
            np.log(sigma_energy),
            np.log(np.maximum(sigma_total, np.finfo(float).tiny)),
            kind="linear",
            bounds_error=False,
            fill_value=(
                float(np.log(max(sigma_total[0], np.finfo(float).tiny))),
                float(np.log(max(sigma_total[-1], np.finfo(float).tiny))),
            ),
        )
        r_interp = interp1d(
            np.log(angle_energy),
            r_theta,
            kind="linear",
            bounds_error=False,
            fill_value=(float(r_theta[0]), float(r_theta[-1])),
        )

        def evaluate(energy_ev: np.ndarray) -> np.ndarray:
            energy_array = np.asarray(energy_ev, dtype=float)
            sigma_e = np.clip(energy_array, sigma_energy[0], sigma_energy[-1])
            angle_e = np.clip(energy_array, angle_energy[0], angle_energy[-1])
            return np.exp(sigma_interp(np.log(sigma_e))) * r_interp(np.log(angle_e))

        return evaluate

    def lookup_expectation(
        self,
        neutral_temperature_ev: float,
        ion_temperature_ev: float,
        *,
        speed_points: int,
        sigma_extent: float,
        tlmt_el_ev_amu: float | None,
    ) -> dict[str, float]:
        speed = speed_grid(neutral_temperature_ev, speed_points, sigma_extent)
        pdf = maxwell_speed_pdf(speed, neutral_temperature_ev)
        safe_speed = np.maximum(speed, 1.0e-30)
        safe_speed2 = safe_speed**2

        # D + D: m_ref = m_i / 2, so the CDF specific-energy and
        # specific-temperature inputs are half the physical lab values.
        table_energy_ev_amu = 0.25 * M_D_KG * safe_speed2 / EV_TO_J
        table_temp_ev_amu = 0.5 * ion_temperature_ev
        table_temp_for_i01 = table_temp_ev_amu
        if tlmt_el_ev_amu is not None:
            table_temp_for_i01 = max(table_temp_for_i01, tlmt_el_ev_amu)
        transport_energy_axis = self.i_1_0.energy_axis
        low_energy_clip = table_energy_ev_amu < transport_energy_axis[0]
        high_energy_clip = table_energy_ev_amu > transport_energy_axis[-1]

        reaction_rate_cm3_s = self.reaction_rate(
            table_energy_ev_amu, table_temp_ev_amu
        )
        i_1_0_m3_s = (
            self.i_1_0(table_energy_ev_amu, table_temp_for_i01) * self.i_1_0_mult
        )
        i_1_1_m4_s2 = (
            self.i_1_1_up(table_energy_ev_amu, table_temp_for_i01)
            * self.i_1_1_mult
        )
        i_1_2_m5_s3 = (
            self.i_1_2_up2(table_energy_ev_amu, table_temp_ev_amu)
            * self.i_1_2_mult
        )

        va2 = 2.0 * ion_temperature_ev * EV_TO_J / M_D_KG
        ame = 0.5
        sp0 = i_1_1_m4_s2 / safe_speed - 0.5 * va2 / safe_speed2 * i_1_0_m3_s
        dtr_eng_j_m3_s = -ame * (
            M_D_KG * sp0 * safe_speed2 - 0.5 * M_D_KG * i_1_2_m5_s3
        )
        dtr_eng_j_m3_s[speed < 1.0e-10] = 0.0

        gain_j_m3_s = float(integrate.simpson(pdf * dtr_eng_j_m3_s, x=speed))
        rate_cm3_s = float(integrate.simpson(pdf * reaction_rate_cm3_s, x=speed))
        gain_ev_cm3_s = gain_j_m3_s / EV_TO_J * 1.0e6
        mean_delta_e_ev = gain_ev_cm3_s / rate_cm3_s if rate_cm3_s > 0.0 else 0.0
        return {
            "gain_ev_cm3_s": gain_ev_cm3_s,
            "loss_ev_cm3_s": -gain_ev_cm3_s,
            "reaction_rate_cm3_s": rate_cm3_s,
            "mean_delta_e_ev": mean_delta_e_ev,
            "i01_table_temperature_ev_amu": table_temp_for_i01,
            "table_temperature_was_clamped": table_temp_for_i01 != table_temp_ev_amu,
            "energy_low_clip_probability": float(
                integrate.simpson(pdf * low_energy_clip.astype(float), x=speed)
            ),
            "energy_high_clip_probability": float(
                integrate.simpson(pdf * high_energy_clip.astype(float), x=speed)
            ),
        }

    def angle_sigma_mt_expectation(
        self,
        neutral_temperature_ev: float,
        ion_temperature_ev: float,
        *,
        speed_points: int,
        sigma_extent: float,
    ) -> float:
        temperature_sum = neutral_temperature_ev + ion_temperature_ev
        if temperature_sum <= 0.0:
            return 0.0

        temp_difference_factor = (neutral_temperature_ev - ion_temperature_ev) / (
            2.0 * temperature_sum
        )
        if temp_difference_factor == 0.0:
            return 0.0

        relative_sigma = np.sqrt(temperature_sum * EV_TO_J / M_D_KG)
        speed = np.linspace(0.0, sigma_extent * relative_sigma, speed_points)
        pdf = (
            np.sqrt(2.0 / np.pi)
            * speed**2
            / relative_sigma**3
            * np.exp(-(speed**2) / (2.0 * relative_sigma**2))
        )
        rel_energy_ev = 0.25 * M_D_KG * speed**2 / EV_TO_J
        sigma_mt_cm2 = self.sigma_mt_from_angle(rel_energy_ev)
        moment = float(integrate.simpson(pdf * sigma_mt_cm2 * speed**3, x=speed))
        prefactor = -(M_D_KG / (2.0 * EV_TO_J)) * SPEED_CM_PER_M
        return prefactor * temp_difference_factor * moment


def build_pairs(args: argparse.Namespace, first_cdf: Path | None) -> list[tuple[float, float]]:
    if args.pairs:
        return parse_pair_list(args.pairs)
    if args.use_transport_temp_grid:
        if first_cdf is None:
            raise ValueError("No CDF file is available for --use-transport-temp-grid.")
        cdf = ElasticCdf.from_file(first_cdf)
        _, temp_axis = transport_axes(cdf)
        temperatures = [2.0 * float(value) for value in temp_axis]
    else:
        temperatures = parse_temperature_list(args.temperatures)
    if args.ion_temperature is not None:
        return [(temperature, float(args.ion_temperature)) for temperature in temperatures]
    return [(temperature, temperature) for temperature in temperatures]


def discover_cdfs(base_dir: Path, args: argparse.Namespace) -> list[Path]:
    if args.all_cdfs:
        return sorted(path for path in base_dir.rglob("*.cdf") if path.is_file())
    candidates = args.cdf or list(DEFAULT_CDFS)
    return [Path(candidate) if Path(candidate).is_absolute() else base_dir / candidate for candidate in candidates]


def format_float(value: float) -> str:
    return f"{value:.6e}"


def print_summary(rows: list[dict[str, object]], summary_temperature: float) -> None:
    equal_rows = [
        row
        for row in rows
        if abs(float(row["neutral_temperature_ev"]) - summary_temperature) < 1.0e-12
        and abs(float(row["ion_temperature_ev"]) - summary_temperature) < 1.0e-12
    ]
    if not equal_rows:
        return

    print("")
    print(f"Equal-temperature summary at Tn = Ti = {summary_temperature:g} eV")
    print("  sign: positive means neutral energy gain")
    for row in equal_rows:
        label = str(row["cdf_path"])
        gain = float(row["lookup_gain_ev_cm3_s"])
        no_clamp = float(row["lookup_gain_no_clamp_ev_cm3_s"])
        clamped_reference = float(row["lookup_gain_clamped_reference_ev_cm3_s"])
        mean_delta = float(row["lookup_mean_delta_e_ev"])
        clamp_state = "on" if row["temperature_clamp_enabled"] else "off"
        print(
            f"  {label}: lookup={format_float(gain)} eV cm^3/s "
            f"(clamp={clamp_state}), no_clamp={format_float(no_clamp)}, "
            f"clamped_ref={format_float(clamped_reference)}, "
            f"mean_dE={format_float(mean_delta)} eV"
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cdf",
        action="append",
        help="CDF/CDL file to evaluate. Can be repeated. Defaults to key el_data CDFs.",
    )
    parser.add_argument(
        "--all-cdfs",
        action="store_true",
        help="Evaluate every *.cdf under el_data instead of the default key CDFs.",
    )
    parser.add_argument(
        "--temperatures",
        nargs="*",
        default=[str(value) for value in DEFAULT_TEMPERATURES_EV],
        help="Physical temperatures in eV. Comma-separated values are accepted.",
    )
    parser.add_argument(
        "--pairs",
        nargs="*",
        help="Explicit physical temperature pairs as T_neutral,T_ion in eV.",
    )
    parser.add_argument(
        "--ion-temperature",
        type=float,
        help="Use a fixed physical ion temperature in eV while scanning neutral temperatures.",
    )
    parser.add_argument(
        "--use-transport-temp-grid",
        action="store_true",
        help="Use physical temperatures equal to 2x the first CDF transport temperature grid.",
    )
    parser.add_argument(
        "--tlmt-el",
        type=float,
        default=0.9,
        help="Fortran TLMT_EL clamp value for I_1_0 and I_1_1, in table eV/amu.",
    )
    clamp_group = parser.add_mutually_exclusive_group()
    clamp_group.add_argument(
        "--temperature-clamp",
        dest="temperature_clamp",
        action="store_true",
        help="Enable the Fortran TLMT_EL clamp for the main lookup result.",
    )
    clamp_group.add_argument(
        "--no-temperature-clamp",
        dest="temperature_clamp",
        action="store_false",
        help="Disable the temperature clamp for the main lookup result.",
    )
    parser.add_argument(
        "--speed-points",
        type=int,
        default=20001,
        help="Number of speed-grid points for Simpson integration.",
    )
    parser.add_argument(
        "--sigma-extent",
        type=float,
        default=12.0,
        help="Integrate speed distributions from 0 to this many thermal sigmas.",
    )
    parser.add_argument(
        "--output",
        default="energy_exchange_bias_summary.csv",
        help="CSV output path.",
    )
    parser.add_argument(
        "--summary-temperature",
        type=float,
        default=2.0,
        help="Equal-temperature case to echo to stdout when present.",
    )
    parser.set_defaults(temperature_clamp=False)
    args = parser.parse_args()

    base_dir = Path(__file__).resolve().parent
    cdf_paths = discover_cdfs(base_dir, args)
    if not cdf_paths:
        raise ValueError("No CDF files were selected.")
    missing = [path for path in cdf_paths if not path.exists()]
    if missing:
        missing_list = ", ".join(str(path) for path in missing)
        raise FileNotFoundError(f"Missing CDF file(s): {missing_list}")

    pairs = build_pairs(args, cdf_paths[0])
    if args.tlmt_el <= 0.0:
        raise ValueError("--tlmt-el must be positive. Use --no-temperature-clamp to disable the clamp.")
    tlmt_el_value = float(args.tlmt_el)
    tlmt_el = tlmt_el_value if args.temperature_clamp else None
    output_path = Path(args.output)
    if not output_path.is_absolute():
        output_path = base_dir / output_path

    rows: list[dict[str, object]] = []
    for cdf_path in cdf_paths:
        diagnostic = EnergyExchangeDiagnostic(cdf_path)
        for neutral_temp, ion_temp in pairs:
            lookup = diagnostic.lookup_expectation(
                neutral_temp,
                ion_temp,
                speed_points=args.speed_points,
                sigma_extent=args.sigma_extent,
                tlmt_el_ev_amu=tlmt_el,
            )
            lookup_no_clamp = diagnostic.lookup_expectation(
                neutral_temp,
                ion_temp,
                speed_points=args.speed_points,
                sigma_extent=args.sigma_extent,
                tlmt_el_ev_amu=None,
            )
            lookup_clamped_reference = diagnostic.lookup_expectation(
                neutral_temp,
                ion_temp,
                speed_points=args.speed_points,
                sigma_extent=args.sigma_extent,
                tlmt_el_ev_amu=tlmt_el_value,
            )
            angle_gain = diagnostic.angle_sigma_mt_expectation(
                neutral_temp,
                ion_temp,
                speed_points=args.speed_points,
                sigma_extent=args.sigma_extent,
            )
            rows.append(
                {
                    "cdf_path": str(cdf_path.relative_to(base_dir)),
                    "neutral_temperature_ev": neutral_temp,
                    "ion_temperature_ev": ion_temp,
                    "table_ion_temperature_ev_amu": 0.5 * ion_temp,
                    "temperature_clamp_enabled": args.temperature_clamp,
                    "tlmt_el_ev_amu": tlmt_el_value,
                    "i01_table_temperature_ev_amu": lookup[
                        "i01_table_temperature_ev_amu"
                    ],
                    "table_temperature_was_clamped": lookup[
                        "table_temperature_was_clamped"
                    ],
                    "energy_low_clip_probability": lookup[
                        "energy_low_clip_probability"
                    ],
                    "energy_high_clip_probability": lookup[
                        "energy_high_clip_probability"
                    ],
                    "lookup_gain_ev_cm3_s": lookup["gain_ev_cm3_s"],
                    "lookup_loss_ev_cm3_s": lookup["loss_ev_cm3_s"],
                    "lookup_reaction_rate_cm3_s": lookup["reaction_rate_cm3_s"],
                    "lookup_mean_delta_e_ev": lookup["mean_delta_e_ev"],
                    "lookup_gain_no_clamp_ev_cm3_s": lookup_no_clamp["gain_ev_cm3_s"],
                    "lookup_mean_delta_no_clamp_e_ev": lookup_no_clamp[
                        "mean_delta_e_ev"
                    ],
                    "lookup_gain_clamped_reference_ev_cm3_s": lookup_clamped_reference[
                        "gain_ev_cm3_s"
                    ],
                    "lookup_mean_delta_clamped_reference_e_ev": lookup_clamped_reference[
                        "mean_delta_e_ev"
                    ],
                    "angle_sigma_mt_gain_ev_cm3_s": angle_gain,
                }
            )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {output_path}")
    print(f"Evaluated {len(cdf_paths)} CDF file(s) and {len(pairs)} temperature pair(s).")
    print_summary(rows, args.summary_temperature)


if __name__ == "__main__":
    main()
