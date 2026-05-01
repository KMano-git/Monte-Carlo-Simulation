#!/usr/bin/env python3
"""Plot equal-temperature elastic energy-exchange bias.

The plotted value is the ion-side energy contribution,

    -<sigma v DeltaE_n>  [eV cm^3/s],

where <sigma v DeltaE_n> is computed by check_energy_exchange_bias.py.  Positive
values therefore mean ion heating, and negative values mean ion cooling.
"""

from __future__ import annotations

import argparse
import csv
import os
import tempfile
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "matplotlib"))

import matplotlib.pyplot as plt
import numpy as np

from check_energy_exchange_bias import EV_TO_J, EnergyExchangeDiagnostic


DEFAULT_OUTPUT = "figure/energy_exchange_ion_equal_temperature.png"
DEFAULT_CSV = "figure/energy_exchange_ion_equal_temperature.csv"
DEFAULT_SIMULATION_CSV = "figure/energy_exchange_simulation_converted.csv"
DEFAULT_SIMULATION_RESULTS = "energy_exchange_bias.txt"
DEFAULT_ION_DENSITY_M3 = 1.0e21
DEFAULT_NEUTRAL_DENSITY_M3 = 5.0e19

KRSTIC_PLOT_CDFS = (
    "Krstic/krstic_dd_total_elastic_integral_priority.cdf",
    "Krstic/krstic_dd_pure_dcs_compat.cdf",
)
FOUR_PLOT_CDFS = (
    "Bachmann/bachmann_dd_from_split_tables_compat.cdf",
    "dd_00_elastic_pure_el_angle_fixed.cdf",
    "Krstic/krstic_dd_total_elastic_integral_priority.cdf",
    "Krstic/krstic_dd_pure_dcs_compat.cdf",
)

PRESET_CDFS = {
    "krstic": KRSTIC_PLOT_CDFS,
    "four": FOUR_PLOT_CDFS,
}

LABELS = {
    "Bachmann/bachmann_dd_from_split_tables_compat.cdf": "Bachmann",
    "dd_00_elastic_pure_el_angle_fixed.cdf": "Bachmann-Janev",
    "Krstic/krstic_dd_total_elastic_integral_priority.cdf": "Krstic total",
    "Krstic/krstic_dd_pure_dcs_compat.cdf": "Krstic pure elastic",
}

SIMULATION_LABELS = {
    "krstic_total": "Krstic_total",
    "krstic_pure_el": "Krstic_pure_el",
}


def resolve_cdf_paths(
    base_dir: Path,
    cdf_args: list[str] | None,
    preset: str,
) -> list[Path]:
    candidates = cdf_args or list(PRESET_CDFS[preset])
    paths = [
        Path(candidate) if Path(candidate).is_absolute() else base_dir / candidate
        for candidate in candidates
    ]
    missing = [path for path in paths if not path.exists()]
    if missing:
        missing_list = ", ".join(str(path) for path in missing)
        raise FileNotFoundError(f"Missing CDF file(s): {missing_list}")
    return paths


def output_path(base_dir: Path, value: str) -> Path:
    path = Path(value)
    if not path.is_absolute():
        path = base_dir / path
    return path


def cdf_label(base_dir: Path, path: Path) -> str:
    try:
        relative = str(path.relative_to(base_dir))
    except ValueError:
        relative = str(path)
    return LABELS.get(relative, relative)


def compute_rows(
    base_dir: Path,
    cdf_paths: list[Path],
    temperatures: np.ndarray,
    *,
    speed_points: int,
    sigma_extent: float,
    tlmt_el_ev_amu: float | None,
    tlmt_el_value: float,
    temperature_clamp_enabled: bool,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for cdf_path in cdf_paths:
        diagnostic = EnergyExchangeDiagnostic(cdf_path)
        label = cdf_label(base_dir, cdf_path)
        relative = str(cdf_path.relative_to(base_dir)) if cdf_path.is_relative_to(base_dir) else str(cdf_path)
        for temperature_ev in temperatures:
            lookup = diagnostic.lookup_expectation(
                float(temperature_ev),
                float(temperature_ev),
                speed_points=speed_points,
                sigma_extent=sigma_extent,
                tlmt_el_ev_amu=tlmt_el_ev_amu,
            )
            neutral_gain = lookup["gain_ev_cm3_s"]
            rows.append(
                {
                    "cdf_path": relative,
                    "label": label,
                    "temperature_ev": float(temperature_ev),
                    "neutral_gain_ev_cm3_s": neutral_gain,
                    "ion_energy_contribution_ev_cm3_s": -neutral_gain,
                    "reaction_rate_cm3_s": lookup["reaction_rate_cm3_s"],
                    "neutral_mean_delta_e_ev": lookup["mean_delta_e_ev"],
                    "ion_mean_delta_e_ev": -lookup["mean_delta_e_ev"],
                    "temperature_clamp_enabled": temperature_clamp_enabled,
                    "tlmt_el_ev_amu": tlmt_el_value,
                    "i01_table_temperature_ev_amu": lookup[
                        "i01_table_temperature_ev_amu"
                    ],
                    "table_temperature_was_clamped": lookup[
                        "table_temperature_was_clamped"
                    ],
                }
            )
    return rows


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def wi_to_ev_cm3_s_factor(ion_density_m3: float, neutral_density_m3: float) -> float:
    return ion_density_m3 * neutral_density_m3 * 1.0e-6 * EV_TO_J


def read_simulation_rows(
    path: Path,
    *,
    ion_density_m3: float,
    neutral_density_m3: float,
) -> list[dict[str, object]]:
    if not path.exists():
        return []

    conversion = wi_to_ev_cm3_s_factor(ion_density_m3, neutral_density_m3)
    rows: list[dict[str, object]] = []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, skipinitialspace=True)
        for raw_row in reader:
            cross_section = str(raw_row["cross-section"]).strip()
            label = SIMULATION_LABELS.get(cross_section, cross_section)
            wi_w_m3 = float(raw_row["Wi"])
            stddev_w_m3 = float(raw_row["stddev"])
            rows.append(
                {
                    "cross_section": cross_section,
                    "label": label,
                    "temperature_ev": float(raw_row["temperature(eV)"]),
                    "wi_w_m3": wi_w_m3,
                    "stddev_w_m3": stddev_w_m3,
                    "ion_energy_contribution_ev_cm3_s": wi_w_m3 / conversion,
                    "stddev_ev_cm3_s": stddev_w_m3 / conversion,
                    "ion_density_m3": ion_density_m3,
                    "neutral_density_m3": neutral_density_m3,
                }
            )
    return rows


def choose_linthresh(values: np.ndarray) -> float:
    finite_abs = np.abs(values[np.isfinite(values)])
    finite_abs = finite_abs[finite_abs > 0.0]
    if finite_abs.size == 0:
        return 1.0e-11
    return max(1.0e-11, float(np.percentile(finite_abs, 5)) * 0.2)


def format_power_tick(value: float) -> str:
    if value == 0.0:
        return "0"
    sign = "-" if value < 0.0 else ""
    exponent = int(round(np.log10(abs(value))))
    return rf"${sign}10^{{{exponent}}}$"


def plot_rows(
    path: Path,
    rows: list[dict[str, object]],
    simulation_rows: list[dict[str, object]],
    *,
    temperature_clamp_enabled: bool,
    tlmt_el_value: float,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    labels = list(dict.fromkeys(str(row["label"]) for row in rows))
    all_plot_values = [float(row["ion_energy_contribution_ev_cm3_s"]) for row in rows]
    all_plot_values.extend(
        float(row["ion_energy_contribution_ev_cm3_s"]) for row in simulation_rows
    )
    all_values = np.array(
        all_plot_values,
        dtype=float,
    )
    linthresh = choose_linthresh(all_values)

    fig, ax = plt.subplots(figsize=(9.0, 5.6))
    color_by_label: dict[str, str] = {}
    for label in labels:
        label_rows = [row for row in rows if row["label"] == label]
        x = np.array([float(row["temperature_ev"]) for row in label_rows])
        y = np.array(
            [float(row["ion_energy_contribution_ev_cm3_s"]) for row in label_rows]
        )
        order = np.argsort(x)
        (line,) = ax.plot(x[order], y[order], lw=2.0, label=label)
        color_by_label[label] = line.get_color()

    simulation_labels = list(
        dict.fromkeys(str(row["label"]) for row in simulation_rows)
    )
    marker_by_label = {
        "Krstic_total": "o",
        "Krstic_pure_el": "s",
    }
    for label in simulation_labels:
        label_rows = [row for row in simulation_rows if row["label"] == label]
        x = np.array([float(row["temperature_ev"]) for row in label_rows])
        y = np.array(
            [float(row["ion_energy_contribution_ev_cm3_s"]) for row in label_rows]
        )
        yerr = np.array([float(row["stddev_ev_cm3_s"]) for row in label_rows])
        order = np.argsort(x)
        ax.errorbar(
            x[order],
            y[order],
            yerr=yerr[order],
            fmt=marker_by_label.get(label, "o"),
            ms=4.8,
            capsize=2.5,
            elinewidth=1.1,
            lw=0.0,
            color=color_by_label.get(label, "0.25"),
            markeredgecolor="white",
            markeredgewidth=0.5,
            label=f"{label} simulation",
        )

    ax.axhline(0.0, color="0.2", lw=0.9)
    ax.set_xscale("log")
    ax.set_yscale("symlog", linthresh=linthresh, linscale=0.8)
    ax.set_xlabel(r"Ion / neutral temperature  $T_i = T_n$  [eV]",fontsize=14)
    ax.set_ylabel(r"Ion energy contribution  $-\langle\sigma v \Delta E_n\rangle$  [eV cm$^3$/s]",fontsize=14)
    # clamp_text = (
    #     f"temperature clamp on (TLMT_EL={tlmt_el_value:g} eV/amu)"
    #     if temperature_clamp_enabled
    #     else "temperature clamp off"
    # )
    # title_suffix = ", with simulation Wi" if simulation_rows else ""
    # ax.set_title(f"Elastic energy contribution bias",fontsize=18)
    ax.grid(True, which="both", alpha=0.28)
    y_min, y_max = ax.get_ylim()
    tick_candidates = np.array(
        [-1.0e-7, -1.0e-8, -1.0e-9, -1.0e-10, -1.0e-11, 0.0, 1.0e-11, 1.0e-10, 1.0e-9, 1.0e-8]
    )
    yticks = [tick for tick in tick_candidates if y_min <= tick <= y_max]
    ax.set_yticks(yticks)
    ax.set_yticklabels([format_power_tick(tick) for tick in yticks])
    ax.legend(fontsize=12)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cdf",
        action="append",
        help="CDF/CDL file to evaluate. Can be repeated. Defaults to Krstic total/pure CDFs.",
    )
    parser.add_argument(
        "--preset",
        choices=sorted(PRESET_CDFS),
        default="krstic",
        help="Default CDF set to plot when --cdf is not specified.",
    )
    parser.add_argument(
        "--temp-min",
        type=float,
        default=0.2,
        help="Minimum equal ion/neutral physical temperature in eV.",
    )
    parser.add_argument(
        "--temp-max",
        type=float,
        default=100.0,
        help="Maximum equal ion/neutral physical temperature in eV.",
    )
    parser.add_argument(
        "--temp-points",
        type=int,
        default=120,
        help="Number of log-spaced temperature points.",
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
        help="Enable the Fortran TLMT_EL clamp for the plotted result.",
    )
    clamp_group.add_argument(
        "--no-temperature-clamp",
        dest="temperature_clamp",
        action="store_false",
        help="Disable the temperature clamp for the plotted result.",
    )
    parser.add_argument(
        "--speed-points",
        type=int,
        default=8001,
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
        default=DEFAULT_OUTPUT,
        help="PNG output path.",
    )
    parser.add_argument(
        "--csv",
        default=DEFAULT_CSV,
        help="CSV output path for plotted values.",
    )
    parser.add_argument(
        "--simulation-results",
        default=DEFAULT_SIMULATION_RESULTS,
        help="Simulation Wi/stddev CSV text to overlay when present.",
    )
    parser.add_argument(
        "--simulation-csv",
        default=DEFAULT_SIMULATION_CSV,
        help="CSV output path for converted simulation values.",
    )
    parser.add_argument(
        "--no-simulation-results",
        action="store_true",
        help="Do not overlay simulation results.",
    )
    parser.add_argument(
        "--ion-density",
        type=float,
        default=DEFAULT_ION_DENSITY_M3,
        help="Ion density used to convert Wi [W/m^3] to eV cm^3/s.",
    )
    parser.add_argument(
        "--neutral-density",
        type=float,
        default=DEFAULT_NEUTRAL_DENSITY_M3,
        help="Neutral density used to convert Wi [W/m^3] to eV cm^3/s.",
    )
    parser.set_defaults(temperature_clamp=False)
    args = parser.parse_args()

    if args.temp_min <= 0.0 or args.temp_max <= 0.0:
        raise ValueError("Temperature bounds must be positive for log spacing.")
    if args.temp_min >= args.temp_max:
        raise ValueError("--temp-min must be smaller than --temp-max.")
    if args.temp_points < 2:
        raise ValueError("--temp-points must be at least 2.")
    if args.tlmt_el <= 0.0:
        raise ValueError("--tlmt-el must be positive. Use --no-temperature-clamp to disable the clamp.")
    if args.ion_density <= 0.0 or args.neutral_density <= 0.0:
        raise ValueError("--ion-density and --neutral-density must be positive.")

    base_dir = Path(__file__).resolve().parent
    cdf_paths = resolve_cdf_paths(base_dir, args.cdf, args.preset)
    temperatures = np.geomspace(args.temp_min, args.temp_max, args.temp_points)
    tlmt_el_value = float(args.tlmt_el)
    tlmt_el = tlmt_el_value if args.temperature_clamp else None

    rows = compute_rows(
        base_dir,
        cdf_paths,
        temperatures,
        speed_points=args.speed_points,
        sigma_extent=args.sigma_extent,
        tlmt_el_ev_amu=tlmt_el,
        tlmt_el_value=tlmt_el_value,
        temperature_clamp_enabled=args.temperature_clamp,
    )
    csv_path = output_path(base_dir, args.csv)
    png_path = output_path(base_dir, args.output)
    simulation_rows: list[dict[str, object]] = []
    if not args.no_simulation_results:
        simulation_path = output_path(base_dir, args.simulation_results)
        simulation_rows = read_simulation_rows(
            simulation_path,
            ion_density_m3=float(args.ion_density),
            neutral_density_m3=float(args.neutral_density),
        )
        plotted_labels = {str(row["label"]) for row in rows}
        simulation_rows = [
            row for row in simulation_rows if str(row["label"]) in plotted_labels
        ]

    write_csv(csv_path, rows)
    if simulation_rows:
        simulation_csv_path = output_path(base_dir, args.simulation_csv)
        write_csv(simulation_csv_path, simulation_rows)
    plot_rows(
        png_path,
        rows,
        simulation_rows,
        temperature_clamp_enabled=args.temperature_clamp,
        tlmt_el_value=tlmt_el_value,
    )

    print(f"Wrote {csv_path}")
    if simulation_rows:
        print(f"Wrote {simulation_csv_path}")
    print(f"Wrote {png_path}")
    print(
        f"Plotted {len(cdf_paths)} CDF file(s), {len(temperatures)} equal-temperature points, "
        f"and {len(simulation_rows)} simulation point(s)."
    )


if __name__ == "__main__":
    main()
