#!/usr/bin/env python3
"""
Fortran (3d3v_event-driven) vs Python (energy_histogram.py) の deltaE ヒストグラム比較スクリプト

Fortranの deltaE は中性粒子のエネルギー変化 (E_new - E_old)
Pythonの deltaE はイオンへのエネルギー移行 (E_n - E_i for CX, E_n - E_n_new for EL)
→ 符号が逆なので、Fortranの符号を反転して比較する
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import re
import os
import sys

# ============================================================
# 物理定数
# ============================================================
AMU = 1.66054e-27
M_D = 2.014 * AMU
EV_TO_J = 1.60218e-19
CM2_TO_M2 = 1e-4

N_SAMPLES = 10000000

# ============================================================
# CDF / 断面積 (energy_histogram.py と同じ)
# ============================================================
def load_elastic_cdf(filename):
    with open(filename, 'r') as f:
        content = f.read()
    match = re.search(r'xs_data_tab\s*=([\d\s.,eE+-]+);', content, re.DOTALL)
    data_str = match.group(1).replace('\n', ' ').replace(',', ' ')
    data = np.fromstring(data_str, sep=' ')
    sigma_raw = data[0:101] * CM2_TO_M2
    angle_data = data[10505:10505+12801]
    angle_cdf = angle_data.reshape((251, 51), order='F')
    energy_grid = np.logspace(np.log10(0.001), np.log10(100), 101)
    energy_grid_angle = np.logspace(np.log10(0.001), np.log10(100), 51)
    prob_grid = np.linspace(0.0, 1.0, 251)
    return energy_grid, sigma_raw, prob_grid, energy_grid_angle, angle_cdf

from scipy.interpolate import RegularGridInterpolator, interp1d

def get_elastic_sigma(energy, energy_grid, sigma_elastic):
    log_interp = interp1d(np.log(energy_grid), np.log(sigma_elastic),
        kind='linear', bounds_error=False,
        fill_value=(np.log(sigma_elastic[0]), np.log(sigma_elastic[-1])))
    return np.exp(log_interp(np.log(np.clip(energy, energy_grid[0], energy_grid[-1]))))

def get_cx_sigma(energy):
    E = np.maximum(energy, 0.01)
    sigma = 0.6937e-18 * (1.0 - 0.155 * np.log10(E))**2 / (1.0 + 0.1112e-14 * E**3.3)
    return np.maximum(sigma, 1e-22)

def sample_maxwell_velocity(temperature, mass, n_samples):
    v_th = np.sqrt(temperature * EV_TO_J / mass)
    return np.random.normal(0, v_th, size=(n_samples, 3))

def sample_scattering_angle_vec(energy, prob_grid, energy_grid_angle, angle_cdf, n):
    rand_prob = np.random.uniform(0, 1, n)
    energy = np.atleast_1d(energy)
    if len(energy) == 1: energy = np.full(n, energy[0])
    energy = np.clip(energy, energy_grid_angle[0], energy_grid_angle[-1])
    interp = RegularGridInterpolator((prob_grid, np.log(energy_grid_angle)),
        angle_cdf, method='linear', bounds_error=False, fill_value=None)
    return interp(np.column_stack([rand_prob, np.log(energy)]))

def rotate_vector(v_rel, chi, phi):
    n = v_rel.shape[0]
    v_mag = np.linalg.norm(v_rel, axis=1, keepdims=True)
    v_mag = np.where(v_mag < 1e-30, 1e-30, v_mag)
    e_z = v_rel / v_mag
    e_x = np.zeros_like(e_z)
    z_al = np.abs(e_z[:, 2]) > 0.9
    e_x[z_al, 0] = -e_z[z_al, 1]; e_x[z_al, 1] = e_z[z_al, 0]
    e_x[~z_al, 0] = -e_z[~z_al, 2]; e_x[~z_al, 2] = e_z[~z_al, 0]
    e_x_m = np.linalg.norm(e_x, axis=1, keepdims=True)
    e_x_m = np.where(e_x_m < 1e-30, 1e-30, e_x_m)
    e_x /= e_x_m; e_y = np.cross(e_z, e_x)
    sc, cc = np.sin(chi)[:,None], np.cos(chi)[:,None]
    sp, cp = np.sin(phi)[:,None], np.cos(phi)[:,None]
    return v_mag * (sc*cp*e_x + sc*sp*e_y + cc*e_z)

# ============================================================
# Fortranのdelta_E_hist.csvを読み込み
# ============================================================
def load_fortran_csv(filename):
    """deltaE_hist.csv を読み込んで bin_center, rate_EL, rate_CX を返す"""
    bc, rel, rcx = [], [], []
    with open(filename, 'r') as f:
        lines = f.readlines()
    # skip header (3 lines)
    for line in lines[3:]:
        line = line.strip()
        if not line: continue
        parts = line.split(',')
        if len(parts) >= 3:
            bc.append(float(parts[0]))
            rel.append(float(parts[1]))
            rcx.append(float(parts[2]))
    return np.array(bc), np.array(rel), np.array(rcx)

# ============================================================
# メイン
# ============================================================
def main():
    # --- Fortranの条件を読み取り ---
    T_ion = 2.0     # from input.nml
    T_neutral = 2.0 # T_init=2eV, init_mode=1 (Maxwell)
    n_i = 1.0e21    # from input.nml
    n_n = 5.0e19    # from input.nml, n_init

    print("=" * 60)
    print("Comparison: Fortran (multi-collision) vs Python (1st collision)")
    print(f"  Tn = {T_neutral} eV, Ti = {T_ion} eV (both Maxwell)")
    print("=" * 60)

    # --- Load CDF ---
    cdf_file = "dd_00_elastic_pure_el_angle.cdf"
    energy_grid, sigma_elastic, prob_grid, energy_grid_angle, angle_cdf = load_elastic_cdf(cdf_file)

    # --- Python 計算 (T_neutral=3eV, Maxwell) ---
    print("Computing Python results (Maxwell neutral T=2eV, T_ion=2eV)...")
    np.random.seed(42)
    mu = M_D / 2
    v_n = sample_maxwell_velocity(T_neutral, M_D, N_SAMPLES)
    v_i = sample_maxwell_velocity(T_ion, M_D, N_SAMPLES)
    v_rel = v_n - v_i
    v_rel_mag = np.linalg.norm(v_rel, axis=1)
    E_rel = 0.5 * mu * v_rel_mag**2 / EV_TO_J
    E_n = 0.5 * M_D * np.sum(v_n**2, axis=1) / EV_TO_J
    E_i = 0.5 * M_D * np.sum(v_i**2, axis=1) / EV_TO_J

    s_cx = get_cx_sigma(E_rel)
    s_el = get_elastic_sigma(E_rel, energy_grid, sigma_elastic)

    dE_cx_py = E_n - E_i  # energy to ions (positive = ions gain)
    chi = sample_scattering_angle_vec(E_rel, prob_grid, energy_grid_angle, angle_cdf, N_SAMPLES)
    phi = np.random.uniform(0, 2*np.pi, N_SAMPLES)
    v_rel_new = rotate_vector(v_rel, chi, phi)
    V_cm = 0.5 * (v_n + v_i)
    v_n_new = V_cm + 0.5 * v_rel_new
    E_n_new = 0.5 * M_D * np.sum(v_n_new**2, axis=1) / EV_TO_J
    dE_el_py = E_n - E_n_new  # energy to ions

    w_cx = s_cx * v_rel_mag
    w_el = s_el * v_rel_mag

    mean_dE_cx_py = np.sum(w_cx * dE_cx_py) / np.sum(w_cx)
    mean_dE_el_py = np.sum(w_el * dE_el_py) / np.sum(w_el)

    print(f"  Python CX weighted mean ΔE (to ions) = {mean_dE_cx_py:+.4f} eV")
    print(f"  Python EL weighted mean ΔE (to ions) = {mean_dE_el_py:+.4f} eV")

    # --- Fortran CSV ---
    print("\nLoading Fortran deltaE_hist.csv ...")
    bc_f, rel_f, rcx_f = load_fortran_csv("deltaE_hist_sol_03.csv")
    bin_width_f = bc_f[1] - bc_f[0]

    # Fortran: deltaE = E_new - E_old (中性粒子の変化)
    # イオンへの移行量は -deltaE → bin_center を反転
    bc_f_flipped = -bc_f

    # Fortranの加重平均ΔE
    total_el_f = np.sum(rel_f) * bin_width_f
    total_cx_f = np.sum(rcx_f) * bin_width_f
    # ΔEの加重平均(符号反転後=イオンへの移行量)
    mean_dE_el_f = np.sum(bc_f_flipped * rel_f * bin_width_f) / total_el_f if total_el_f > 0 else 0
    mean_dE_cx_f = np.sum(bc_f_flipped * rcx_f * bin_width_f) / total_cx_f if total_cx_f > 0 else 0
    print(f"  Fortran EL weighted mean ΔE (to ions, sign-flipped) = {mean_dE_el_f:+.4f} eV")
    print(f"  Fortran CX weighted mean ΔE (to ions, sign-flipped) = {mean_dE_cx_f:+.4f} eV")

    # --- ヒストグラム比較プロット ---
    os.makedirs('figure', exist_ok=True)

    # Python: σv重み付きヒストグラムを計算
    hist_range = (-20, 20)
    n_bins = 400

    counts_cx_py, edges = np.histogram(dE_cx_py, bins=n_bins, range=hist_range, weights=w_cx)
    counts_el_py, _ = np.histogram(dE_el_py, bins=n_bins, range=hist_range, weights=w_el)
    bc_py = 0.5 * (edges[:-1] + edges[1:])
    bw_py = edges[1] - edges[0]

    # 正規化: PythonのヒストグラムをFortranと同じ単位 [events / (m^3 * s)] に変換する
    # Pythonでの重みは w = sigma * v_rel。その平均 <sigma v> は sum(w)/N_SAMPLES
    # 物理的な反応率密度 R = n_n * n_i * <sigma v>
    # したがって、各ビンの値 sum(w_bin) に (n_n * n_i / N_SAMPLES) を掛ければよい
    norm_py = n_n * n_i / N_SAMPLES
    
    rate_el_py = counts_el_py * norm_py
    rate_cx_py = counts_cx_py * norm_py

    fig, axes = plt.subplots(2, 1, figsize=(10, 10))

    # --- EL比較 ---
    ax = axes[0]
    ax.step(bc_py, rate_el_py, where='mid', color='blue', linewidth=1.5,
            label=f'Python (1st coll, Maxwell T={T_neutral}eV)\nmean={mean_dE_el_py:+.3f} eV\nTotal Rate={np.sum(rate_el_py):.2e} m⁻³s⁻¹', alpha=0.8)
    # Fortranは符号反転して表示
    ax.step(-bc_f, rel_f, where='mid', color='red', linewidth=1.5, linestyle='--',
            label=f'Fortran (multi-coll, Maxwell T={T_neutral}eV)\nmean={mean_dE_el_f:+.3f} eV\nTotal Rate={np.sum(rel_f):.2e} m⁻³s⁻¹', alpha=0.8)
    ax.axvline(x=0, color='black', linewidth=0.5)
    ax.axvline(x=mean_dE_el_py, color='blue', linewidth=1, linestyle=':')
    ax.axvline(x=mean_dE_el_f, color='red', linewidth=1, linestyle=':')
    ax.set_xlabel('ΔE to ions [eV]', fontsize=13)
    ax.set_ylabel('Reaction Rate Density [m⁻³ s⁻¹]', fontsize=13)
    ax.set_title(f'Elastic Scattering (EL) — Tn={T_neutral}eV, Ti={T_ion}eV', fontsize=14, fontweight='bold')
    ax.set_yscale('log')
    # 自動スケールに任せるが、下限はある程度切る
    max_val = max(np.max(rate_el_py), np.max(rel_f))
    ax.set_ylim(max_val * 1e-5, max_val * 2)
    ax.set_xlim(-15, 15)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)

    # --- CX比較 ---
    ax = axes[1]
    ax.step(bc_py, rate_cx_py, where='mid', color='blue', linewidth=1.5,
            label=f'Python (1st coll, Maxwell T={T_neutral}eV)\nmean={mean_dE_cx_py:+.3f} eV\nTotal Rate={np.sum(rate_cx_py):.2e} m⁻³s⁻¹', alpha=0.8)
    ax.step(-bc_f, rcx_f, where='mid', color='red', linewidth=1.5, linestyle='--',
            label=f'Fortran (multi-coll, Maxwell T={T_neutral}eV)\nmean={mean_dE_cx_f:+.3f} eV\nTotal Rate={np.sum(rcx_f):.2e} m⁻³s⁻¹', alpha=0.8)
    ax.axvline(x=0, color='black', linewidth=0.5)
    ax.axvline(x=mean_dE_cx_py, color='blue', linewidth=1, linestyle=':')
    ax.axvline(x=mean_dE_cx_f, color='red', linewidth=1, linestyle=':')
    ax.set_xlabel('ΔE to ions [eV]', fontsize=13)
    ax.set_ylabel('Reaction Rate Density [m⁻³ s⁻¹]', fontsize=13)
    ax.set_title(f'Charge Exchange (CX) — Tn={T_neutral}eV, Ti={T_ion}eV', fontsize=14, fontweight='bold')
    ax.set_yscale('log')
    max_val_cx = max(np.max(rate_cx_py), np.max(rcx_f))
    ax.set_ylim(1e22, max_val_cx * 2)
    ax.set_xlim(-15, 15)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    outfile = 'figure/deltaE_comparison.png'
    plt.savefig(outfile, dpi=150, bbox_inches='tight')
    print(f"\nComparison plot saved to: {outfile}")

if __name__ == "__main__":
    main()
