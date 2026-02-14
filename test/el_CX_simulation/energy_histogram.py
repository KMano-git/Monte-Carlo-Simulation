#!/usr/bin/env python3
"""
Energy Transfer Histogram Visualization
固定温度でのエネルギー受け渡し量のヒストグラム表示

Usage:
    python energy_histogram.py [Tn] [Ti]
    
    Tn: 中性粒子温度 [eV] (デフォルト: 1.0)
    Ti: イオン温度 [eV] (デフォルト: 5.0)
"""

import numpy as np
import re
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator, interp1d
import sys
import os

# =============================================================================
# 物理定数（SI単位系）
# =============================================================================
AMU = 1.66054e-27       # kg/amu
M_D = 2.014 * AMU       # 重水素質量 [kg]
EV_TO_J = 1.60218e-19   # J/eV
CM2_TO_M2 = 1e-4        # cm² → m²

# シミュレーションパラメータ
N_SAMPLES = 100000      # サンプル数（ヒストグラム用に多めに取る）


# =============================================================================
# CDF ファイル読み込み（メインスクリプトと同じ）
# =============================================================================
def load_elastic_cdf(filename):
    """dd_00_elastic.cdf からデータを読み込む"""
    with open(filename, 'r') as f:
        content = f.read()
    
    match = re.search(r'xs_data_tab\s*=\s*([\d\s.,eE+-]+);', content, re.DOTALL)
    if not match:
        raise ValueError("xs_data_tab not found in CDF file")
    
    data_str = match.group(1)
    data_str = data_str.replace('\n', ' ').replace(',', ' ')
    data = np.fromstring(data_str, sep=' ')
    
    sigma_raw = data[0:101] * CM2_TO_M2
    angle_data = data[10505:10505+12801]
    angle_cdf = angle_data.reshape((251, 51), order='F')
    
    energy_grid = np.logspace(np.log10(0.001), np.log10(100), 101)
    energy_grid_angle = np.logspace(np.log10(0.001), np.log10(100), 51)
    prob_grid = np.linspace(0.0, 1.0, 251)
    
    return energy_grid, sigma_raw, prob_grid, energy_grid_angle, angle_cdf


def get_elastic_sigma(energy, energy_grid, sigma_elastic):
    """弾性散乱断面積を補間で取得"""
    log_interp = interp1d(
        np.log(energy_grid), np.log(sigma_elastic),
        kind='linear', bounds_error=False,
        fill_value=(np.log(sigma_elastic[0]), np.log(sigma_elastic[-1]))
    )
    return np.exp(log_interp(np.log(np.clip(energy, energy_grid[0], energy_grid[-1]))))


def get_cx_sigma(energy):
    """荷電交換断面積（理論式）"""
    E = np.maximum(energy, 0.01)
    sigma = 0.6937e-18 * (1.0 - 0.155 * np.log10(E))**2 / (1.0 + 0.1112e-14 * E**3.3)
    return np.maximum(sigma, 1e-22)


def sample_maxwell_velocity(temperature, mass, n_samples):
    """マクスウェル分布から速度をサンプリング"""
    v_th = np.sqrt(temperature * EV_TO_J / mass)
    velocity = np.random.normal(0, v_th, size=(n_samples, 3))
    return velocity


def sample_scattering_angle_vectorized(energy, prob_grid, energy_grid_angle, angle_cdf, n_samples):
    """CDFから散乱角をサンプリング"""
    rand_prob = np.random.uniform(0, 1, n_samples)
    energy = np.atleast_1d(energy)
    if len(energy) == 1:
        energy = np.full(n_samples, energy[0])
    energy = np.clip(energy, energy_grid_angle[0], energy_grid_angle[-1])
    
    interpolator = RegularGridInterpolator(
        (prob_grid, np.log(energy_grid_angle)),
        angle_cdf,
        method='linear',
        bounds_error=False,
        fill_value=None
    )
    points = np.column_stack([rand_prob, np.log(energy)])
    chi = interpolator(points)
    return chi


def rotate_vector(v_rel, chi, phi):
    """相対速度ベクトルを散乱角で回転"""
    n = v_rel.shape[0]
    v_rel_mag = np.linalg.norm(v_rel, axis=1, keepdims=True)
    v_rel_mag = np.where(v_rel_mag < 1e-30, 1e-30, v_rel_mag)
    e_z = v_rel / v_rel_mag
    e_x = np.zeros_like(e_z)
    z_aligned = np.abs(e_z[:, 2]) > 0.9
    
    e_x[z_aligned, 0] = -e_z[z_aligned, 1]
    e_x[z_aligned, 1] = e_z[z_aligned, 0]
    e_x[z_aligned, 2] = 0
    e_x[~z_aligned, 0] = -e_z[~z_aligned, 2]
    e_x[~z_aligned, 1] = 0
    e_x[~z_aligned, 2] = e_z[~z_aligned, 0]
    
    e_x_mag = np.linalg.norm(e_x, axis=1, keepdims=True)
    e_x_mag = np.where(e_x_mag < 1e-30, 1e-30, e_x_mag)
    e_x = e_x / e_x_mag
    e_y = np.cross(e_z, e_x)
    
    sin_chi = np.sin(chi)[:, np.newaxis]
    cos_chi = np.cos(chi)[:, np.newaxis]
    cos_phi = np.cos(phi)[:, np.newaxis]
    sin_phi = np.sin(phi)[:, np.newaxis]
    
    v_rel_new = v_rel_mag * (
        sin_chi * cos_phi * e_x +
        sin_chi * sin_phi * e_y +
        cos_chi * e_z
    )
    return v_rel_new


def compute_energy_transfers(T_neutral, T_ion, n_samples,
                              energy_grid, sigma_elastic,
                              prob_grid, energy_grid_angle, angle_cdf):
    """
    指定温度でのエネルギー受け渡し量を計算
    
    Returns:
        delta_E_cx: CXによるエネルギー変化 [eV] (n_samples,)
        delta_E_el: 弾性散乱によるエネルギー変化 [eV] (n_samples,)
        sigma_cx: CX断面積 [m²]
        sigma_el: 弾性散乱断面積 [m²]
        v_rel_mag: 相対速度 [m/s]
    """
    mu = M_D / 2
    
    # 速度サンプリング
    v_n = sample_maxwell_velocity(T_neutral, M_D, n_samples)
    v_i = sample_maxwell_velocity(T_ion, M_D, n_samples)
    
    # 相対速度
    v_rel = v_n - v_i
    v_rel_mag = np.linalg.norm(v_rel, axis=1)
    
    # 衝突エネルギー [eV]
    E_rel = 0.5 * mu * v_rel_mag**2 / EV_TO_J
    
    # 初期エネルギー [eV]
    E_n = 0.5 * M_D * np.sum(v_n**2, axis=1) / EV_TO_J
    E_i = 0.5 * M_D * np.sum(v_i**2, axis=1) / EV_TO_J
    
    # ===== 荷電交換 (CX) =====
    sigma_cx = get_cx_sigma(E_rel)
    delta_E_cx = E_i - E_n  # 速度が入れ替わる
    
    # ===== 弾性散乱 (EL) =====
    sigma_el = get_elastic_sigma(E_rel, energy_grid, sigma_elastic)
    
    chi = sample_scattering_angle_vectorized(
        E_rel, prob_grid, energy_grid_angle, angle_cdf, n_samples
    )
    phi = np.random.uniform(0, 2*np.pi, n_samples)
    
    v_rel_new = rotate_vector(v_rel, chi, phi)
    V_cm = 0.5 * (v_n + v_i)
    v_n_new = V_cm + 0.5 * v_rel_new
    
    E_n_new = 0.5 * M_D * np.sum(v_n_new**2, axis=1) / EV_TO_J
    delta_E_el = E_n_new - E_n
    
    return delta_E_cx, delta_E_el, sigma_cx, sigma_el, v_rel_mag


def main():
    # コマンドライン引数から温度を取得
    T_neutral = 1.0  # デフォルト
    T_ion = 5.0      # デフォルト
    
    if len(sys.argv) >= 2:
        try:
            T_neutral = float(sys.argv[1])
        except ValueError:
            print(f"Warning: Invalid Tn value '{sys.argv[1]}', using default {T_neutral}")
    
    if len(sys.argv) >= 3:
        try:
            T_ion = float(sys.argv[2])
        except ValueError:
            print(f"Warning: Invalid Ti value '{sys.argv[2]}', using default {T_ion}")
    
    print("=" * 60)
    print("Energy Transfer Histogram")
    print("=" * 60)
    print(f"  Neutral temperature Tn = {T_neutral:.2f} eV")
    print(f"  Ion temperature Ti = {T_ion:.2f} eV")
    print(f"  Number of samples = {N_SAMPLES}")
    
    # CDFファイルを読み込む
    print("\nLoading CDF file...")
    cdf_file = "dd_00_elastic.cdf"
    energy_grid, sigma_elastic, prob_grid, energy_grid_angle, angle_cdf = \
        load_elastic_cdf(cdf_file)
    
    # エネルギー受け渡しを計算
    print("Computing energy transfers...")
    delta_E_cx, delta_E_el, sigma_cx, sigma_el, v_rel_mag = compute_energy_transfers(
        T_neutral, T_ion, N_SAMPLES,
        energy_grid, sigma_elastic,
        prob_grid, energy_grid_angle, angle_cdf
    )
    
    # 反応率で重み付けした平均
    # <sigma * v * ΔE> / <sigma * v>
    weight_cx = sigma_cx * v_rel_mag
    weight_el = sigma_el * v_rel_mag
    
    total_rate_cx = np.sum(weight_cx)
    total_rate_el = np.sum(weight_el)
    
    mean_delta_E_cx = np.sum(weight_cx * delta_E_cx) / total_rate_cx
    mean_delta_E_el = np.sum(weight_el * delta_E_el) / total_rate_el
    
    # 単純平均も計算
    simple_mean_cx = np.mean(delta_E_cx)
    simple_mean_el = np.mean(delta_E_el)
    
    print(f"\n===== Results =====")
    print(f"Charge Exchange (CX):")
    print(f"  Simple mean ΔE = {simple_mean_cx:+.4f} eV")
    print(f"  Weighted mean ΔE = {mean_delta_E_cx:+.4f} eV")
    print(f"  Total Rate (arb.) = {total_rate_cx:.4e}")
    print(f"Elastic Scattering (EL):")
    print(f"  Simple mean ΔE = {simple_mean_el:+.4f} eV")
    print(f"  Weighted mean ΔE = {mean_delta_E_el:+.4f} eV")
    print(f"  Total Rate (arb.) = {total_rate_el:.4e}")
    
    # プロット切り替え: True=1つの図に統合, False=2つのサブプロット
    COMBINED_PLOT = True
    
    # ヒストグラムの範囲を決定
    all_data = np.concatenate([delta_E_cx, delta_E_el])
    hist_range = (np.percentile(all_data, 1), np.percentile(all_data, 99))
    bins = 80
    
    if COMBINED_PLOT:
        # ===== 統合プロット（1つの図に両方） =====
        fig, ax = plt.subplots(1, 1, figsize=(10, 7))
        
        # CXヒストグラム
        counts_cx, bin_edges_cx, _ = ax.hist(
            delta_E_cx, bins=bins, range=hist_range, weights=weight_cx,
            color='blue', alpha=0.6, edgecolor='darkblue', linewidth=0.5,
            label='CX rate'
        )
        
        # ELヒストグラム
        counts_el, bin_edges_el, _ = ax.hist(
            delta_E_el, bins=bins, range=hist_range, weights=weight_el,
            color='red', alpha=0.6, edgecolor='darkred', linewidth=0.5,
            label='EL rate'
        )
        
        # 平均値の線
        ax.axvline(x=mean_delta_E_cx, color='blue', linestyle='--', linewidth=2,
                   label=f'CX Mean = {mean_delta_E_cx:+.3f} eV')
        ax.axvline(x=mean_delta_E_el, color='red', linestyle='--', linewidth=2,
                   label=f'EL Mean = {mean_delta_E_el:+.3f} eV')
        ax.axvline(x=0, color='black', linestyle='-', linewidth=1, label='ΔE = 0')
        
        ax.set_xlabel('Energy Transfer ΔE [eV]', fontsize=12)
        ax.set_ylabel('Reaction Rate (σv) [arb. units]', fontsize=12)
        ax.set_yscale('log')
        ax.set_title(f'Energy Transfer Comparison: CX vs EL\nTn = {T_neutral:.2f} eV, Ti = {T_ion:.2f} eV', 
                     fontsize=13, fontweight='bold')
        ax.legend(loc='upper right', fontsize=10)
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        # 保存
        output_file = f'figure/energy_histogram_combined_Tn{T_neutral:.1f}_Ti{T_ion:.1f}.png'
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"\nPlot saved to: {output_file}")
        
    else:
        # ===== 個別プロット（2つのサブプロット） =====
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        
        # ===== CX ヒストグラム =====
        ax1 = axes[0]
        # weightsを指定して断面積×速度を反映
        counts_cx, bin_edges_cx, _ = ax1.hist(
            delta_E_cx, bins=bins, range=hist_range, weights=weight_cx,
            color='blue', alpha=0.7, edgecolor='darkblue', linewidth=0.5,
            label='CX rate'
        )
        
        ax1.axvline(x=mean_delta_E_cx, color='red', linestyle='--', linewidth=2,
                    label=f'Weighted Mean = {mean_delta_E_cx:+.3f} eV')
        ax1.axvline(x=0, color='black', linestyle='-', linewidth=1)
        
        ax1.set_xlabel('Energy Transfer ΔE [eV]', fontsize=12)
        ax1.set_ylabel('Reaction Rate (σv) [arb. units]', fontsize=12)
        ax1.set_yscale('log')
        ax1.set_title(f'Charge Exchange (CX)\nTn = {T_neutral:.2f} eV, Ti = {T_ion:.2f} eV', fontsize=12)
        ax1.legend(loc='upper right', fontsize=10)
        ax1.grid(True, alpha=0.3)
        
        # 注釈
        if mean_delta_E_cx > 0:
            ax1.annotate(f'Heating\n(neutral gains energy)', 
                        xy=(0.95, 0.85), xycoords='axes fraction',
                        fontsize=9, ha='right', color='darkgreen')
        else:
            ax1.annotate(f'Cooling\n(neutral loses energy)', 
                        xy=(0.95, 0.85), xycoords='axes fraction',
                        fontsize=9, ha='right', color='darkred')
        
        # ===== EL ヒストグラム =====
        ax2 = axes[1]
        # weightsを指定して断面積×速度を反映
        counts_el, bin_edges_el, _ = ax2.hist(
            delta_E_el, bins=bins, range=hist_range, weights=weight_el,
            color='red', alpha=0.7, edgecolor='darkred', linewidth=0.5,
            label='EL rate'
        )
        
        ax2.axvline(x=mean_delta_E_el, color='blue', linestyle='--', linewidth=2,
                    label=f'Weighted Mean = {mean_delta_E_el:+.3f} eV')
        ax2.axvline(x=0, color='black', linestyle='-', linewidth=1)
        
        ax2.set_xlabel('Energy Transfer ΔE [eV]', fontsize=12)
        ax2.set_ylabel('Reaction Rate (σv) [arb. units]', fontsize=12)
        ax2.set_yscale('log')
        ax2.set_title(f'Elastic Scattering (EL)\nTn = {T_neutral:.2f} eV, Ti = {T_ion:.2f} eV', fontsize=12)
        ax2.legend(loc='upper right', fontsize=10)
        ax2.grid(True, alpha=0.3)
        
        if mean_delta_E_el > 0:
            ax2.annotate(f'Heating\n(neutral gains energy)', 
                        xy=(0.95, 0.85), xycoords='axes fraction',
                        fontsize=9, ha='right', color='darkgreen')
        else:
            ax2.annotate(f'Cooling\n(neutral loses energy)', 
                        xy=(0.95, 0.85), xycoords='axes fraction',
                        fontsize=9, ha='right', color='darkred')
        
        plt.tight_layout()
        
        # 保存
        output_file = f'figure/energy_histogram_Tn{T_neutral:.1f}_Ti{T_ion:.1f}.png'
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"\nPlot saved to: {output_file}")
    
    plt.show()


if __name__ == "__main__":
    main()
