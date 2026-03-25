#!/usr/bin/env python3
"""Plot 1D3V history/profile CSV outputs."""

from __future__ import annotations

import csv
import os
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RUN = ROOT / "run"
FIG = RUN / "figure"
MPLCONFIG = Path(tempfile.gettempdir()) / "monte_simulations_ai_1d3v_mplconfig"

os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIG))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def read_csv_dict(path: Path) -> dict[str, np.ndarray]:
    with path.open("r", newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    data: dict[str, list[float]] = {name: [] for name in reader.fieldnames or []}
    for row in rows:
        for key, value in row.items():
            data[key].append(float(value))
    return {key: np.asarray(vals, dtype=float) for key, vals in data.items()}


def main() -> None:
    history_path = RUN / "ntscrg.csv"
    profile_path = RUN / "profile_x.csv"
    if not history_path.exists() or not profile_path.exists():
        raise SystemExit("run/ntscrg.csv and run/profile_x.csv are required")

    FIG.mkdir(parents=True, exist_ok=True)

    hist = read_csv_dict(history_path)
    prof = read_csv_dict(profile_path)
    has_lookup = "TLLu_Q_total[W/m3]" in hist and "TLLu_Q_total[W/m3]" in prof

    plt.rcParams["figure.dpi"] = 150
    plt.rcParams["font.size"] = 11

    time_us = hist["time[s]"] * 1.0e6

    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    axes[0].plot(time_us, hist["n_alive"], label="n_alive", color="#0B6E4F", linewidth=1.8)
    axes[0].plot(time_us, hist["weight_sum"], label="weight_sum", color="#C84C09", linewidth=1.5)
    axes[0].set_ylabel("Count / Weight")
    axes[0].set_title("1D3V Time History")
    axes[0].grid(True, alpha=0.3)
    axes[0].legend()

    axes[1].plot(time_us, hist["CL_Q_total[W/m3]"], label="CL total", color="#004E98", linewidth=1.6)
    axes[1].plot(time_us, hist["TL_Q_total[W/m3]"], label="TL total", color="#E63946", linewidth=1.6)
    if has_lookup:
        axes[1].plot(time_us, hist["TLLu_Q_total[W/m3]"], label="TLLu total", color="#6A4C93", linewidth=1.4)
    axes[1].set_xlabel("Time [us]")
    axes[1].set_ylabel("Q [W/m3]")
    axes[1].grid(True, alpha=0.3)
    axes[1].legend()

    fig.tight_layout()
    fig.savefig(FIG / "1d3v_summary.png", bbox_inches="tight")
    plt.close(fig)

    x_cm = prof["x_center[m]"] * 100.0

    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    axes[0].plot(x_cm, prof["CL_Q_total[W/m3]"], label="CL total", color="#004E98", linewidth=1.6)
    axes[0].plot(x_cm, prof["TL_Q_total[W/m3]"], label="TL total", color="#E63946", linewidth=1.6)
    if has_lookup:
        axes[0].plot(x_cm, prof["TLLu_Q_total[W/m3]"], label="TLLu total", color="#6A4C93", linewidth=1.4)
    axes[0].set_ylabel("Q [W/m3]")
    axes[0].set_title("1D3V X Profiles")
    axes[0].grid(True, alpha=0.3)
    axes[0].legend()

    axes[1].plot(x_cm, prof["neutral_density[m-3]"], color="#0B6E4F", linewidth=1.8)
    axes[1].set_xlabel("x [cm]")
    axes[1].set_ylabel("n_n [m-3]")
    axes[1].grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(FIG / "1d3v_profiles.png", bbox_inches="tight")
    plt.close(fig)

    print(f"Saved: {FIG / '1d3v_summary.png'}")
    print(f"Saved: {FIG / '1d3v_profiles.png'}")


if __name__ == "__main__":
    main()
