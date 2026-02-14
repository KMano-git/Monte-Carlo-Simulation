#!/usr/bin/env python3
"""
Monte Carlo Simulation for Neutral Particle-Ion Energy Transport
重水素原子(D⁰)と重水素イオン(D⁺)の衝突によるエネルギー輸送率計算

Based on DEGAS2 logic with improved CX cross-section formula
"""

import numpy as np
import re
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d, RegularGridInterpolator

# =============================================================================
# 物理定数（SI単位系）
# =============================================================================
AMU = 1.66054e-27       # kg/amu
M_D = 2.014 * AMU       # 重水素質量 [kg]
EV_TO_J = 1.60218e-19   # J/eV
CM2_TO_M2 = 1e-4        # cm² → m²

# シミュレーションパラメータ
N_SAMPLES = 50000       # 各温度点での試行回数
T_ION = 5.0             # イオン温度 [eV]（固定）


# =============================================================================
# CDF ファイル読み込み
# =============================================================================
def load_elastic_cdf(filename):
    """
    dd_00_elastic.cdf からデータを読み込む
    
    Returns:
        energy_grid: エネルギーグリッド [eV] (101点)
        sigma_elastic: 弾性散乱断面積 [m²] (101点)
        prob_grid: 確率グリッド (251点)
        energy_grid_angle: 散乱角用エネルギーグリッド [eV] (51点)
        angle_cdf: 散乱角CDF配列 (251, 51) - angle_cdf[prob_idx, energy_idx]
    """
    with open(filename, 'r') as f:
        content = f.read()
    
    # xs_data_tab の数値データを抽出
    match = re.search(r'xs_data_tab\s*=\s*([\d\s.,eE+-]+);', content, re.DOTALL)
    if not match:
        raise ValueError("xs_data_tab not found in CDF file")
    
    # 数値をパース
    data_str = match.group(1)
    data_str = data_str.replace('\n', ' ').replace(',', ' ')
    data = np.fromstring(data_str, sep=' ')
    
    # 断面積データ（インデックス 0-100、101点）
    # 単位: cm² → m²
    sigma_raw = data[0:101] * CM2_TO_M2
    
    # 散乱角CDFデータ（インデックス 10505-23305、12801点 = 251×51）
    # Column-major (Fortran) 形式: (確率251 × エネルギー51)
    angle_data = data[10505:10505+12801]
    angle_cdf = angle_data.reshape((251, 51), order='F')
    
    # エネルギーグリッドの構築
    # 仕様書より: 0.001 eV 〜 100 eV（対数等間隔）
    # 断面積用: 101点
    energy_grid = np.logspace(np.log10(0.001), np.log10(100), 101)
    
    # 散乱角用: 51点
    energy_grid_angle = np.logspace(np.log10(0.001), np.log10(100), 51)
    
    # 確率グリッド: 0 → 1 (251点)
    prob_grid = np.linspace(0.0, 1.0, 251)
    
    return energy_grid, sigma_raw, prob_grid, energy_grid_angle, angle_cdf


def get_elastic_sigma(energy, energy_grid, sigma_elastic):
    """
    エネルギーに対応する弾性散乱断面積を補間で取得
    
    Args:
        energy: 衝突エネルギー [eV]
        energy_grid: エネルギーグリッド [eV]
        sigma_elastic: 断面積データ [m²]
    
    Returns:
        断面積 [m²]
    """
    # 対数補間
    log_interp = interp1d(
        np.log(energy_grid), np.log(sigma_elastic),
        kind='linear', bounds_error=False,
        fill_value=(np.log(sigma_elastic[0]), np.log(sigma_elastic[-1]))
    )
    return np.exp(log_interp(np.log(np.clip(energy, energy_grid[0], energy_grid[-1]))))


# =============================================================================
# 荷電交換断面積（理論式）
# =============================================================================
def get_cx_sigma(energy):
    """
    荷電交換断面積を理論式から計算（Janev/Rabinovitch近似）
    
    Args:
        energy: 衝突エネルギー [eV]
    
    Returns:
        断面積 [m²]
    """
    # 簡易フィッティング式
    # 低エネルギーでの発散を防ぐためにクリップ
    E = np.maximum(energy, 0.01)
    sigma = 3.0e-19 * (1.0 - 0.05 * np.log10(E))**2
    return np.maximum(sigma, 1e-22)  # 最小値を設定


# =============================================================================
# 速度サンプリング
# =============================================================================
def sample_maxwell_velocity(temperature, mass, n_samples):
    """
    マクスウェル分布から速度ベクトルをサンプリング
    
    Args:
        temperature: 温度 [eV]
        mass: 質量 [kg]
        n_samples: サンプル数
    
    Returns:
        velocity: 速度ベクトル配列 (n_samples, 3) [m/s]
    """
    # 熱速度 v_th = sqrt(kT/m)、ここで kT は eV 単位
    v_th = np.sqrt(temperature * EV_TO_J / mass)
    
    # 各成分はガウス分布 N(0, v_th)
    velocity = np.random.normal(0, v_th, size=(n_samples, 3))
    
    return velocity


# =============================================================================
# 散乱角サンプリング
# =============================================================================
def sample_scattering_angle(energy, prob_grid, energy_grid_angle, angle_cdf, n_samples):
    """
    CDFから散乱角をサンプリング
    
    Args:
        energy: 衝突エネルギー [eV] (スカラーまたは配列)
        prob_grid: 確率グリッド (251点)
        energy_grid_angle: エネルギーグリッド (51点)
        angle_cdf: 散乱角CDF (251, 51)
        n_samples: サンプル数
    
    Returns:
        chi: 散乱角 [rad] (n_samples,)
    """
    # ランダム確率を生成
    rand_prob = np.random.uniform(0, 1, n_samples)
    
    # エネルギーをクリップ
    energy = np.atleast_1d(energy)
    if len(energy) == 1:
        energy = np.full(n_samples, energy[0])
    energy = np.clip(energy, energy_grid_angle[0], energy_grid_angle[-1])
    
    # 各サンプルについて散乱角を補間
    chi = np.zeros(n_samples)
    
    # エネルギーインデックスを計算（対数補間）
    log_E = np.log(energy)
    log_E_grid = np.log(energy_grid_angle)
    
    for i in range(n_samples):
        # エネルギー補間インデックス
        e_idx = np.searchsorted(log_E_grid, log_E[i]) - 1
        e_idx = np.clip(e_idx, 0, len(energy_grid_angle) - 2)
        
        # エネルギー方向の線形補間係数
        t = (log_E[i] - log_E_grid[e_idx]) / (log_E_grid[e_idx + 1] - log_E_grid[e_idx])
        t = np.clip(t, 0, 1)
        
        # 2つのエネルギー点でのCDFを補間で結合
        cdf_low = angle_cdf[:, e_idx]
        cdf_high = angle_cdf[:, e_idx + 1]
        cdf_interp = (1 - t) * cdf_low + t * cdf_high
        
        # CDFから散乱角を逆変換
        # 確率 0→1 に対して散乱角 π→0 (cos χ = -1 → 1)
        # CDFの値がラジアンでの散乱角を直接表している
        p_idx = np.searchsorted(prob_grid, rand_prob[i]) - 1
        p_idx = np.clip(p_idx, 0, len(prob_grid) - 2)
        
        # 確率方向の線形補間
        s = (rand_prob[i] - prob_grid[p_idx]) / (prob_grid[p_idx + 1] - prob_grid[p_idx])
        s = np.clip(s, 0, 1)
        
        chi[i] = (1 - s) * cdf_interp[p_idx] + s * cdf_interp[p_idx + 1]
    
    return chi


def sample_scattering_angle_vectorized(energy, prob_grid, energy_grid_angle, angle_cdf, n_samples):
    """
    CDFから散乱角をサンプリング（ベクトル化版、高速）
    
    散乱角CDFのデータ:
    - 確率軸 P: 0.0 → 1.0
    - 対応する散乱角 χ: π → 0（ラジアン）
    """
    # ランダム確率を生成
    rand_prob = np.random.uniform(0, 1, n_samples)
    
    # エネルギーをクリップ
    energy = np.atleast_1d(energy)
    if len(energy) == 1:
        energy = np.full(n_samples, energy[0])
    energy = np.clip(energy, energy_grid_angle[0], energy_grid_angle[-1])
    
    # 2D補間器を構築
    # RegularGridInterpolatorは (prob, energy) の順で点を指定
    interpolator = RegularGridInterpolator(
        (prob_grid, np.log(energy_grid_angle)),
        angle_cdf,
        method='linear',
        bounds_error=False,
        fill_value=None
    )
    
    # 補間点を構築
    points = np.column_stack([rand_prob, np.log(energy)])
    chi = interpolator(points)
    
    return chi


# =============================================================================
# 座標変換（ロドリゲスの回転公式）
# =============================================================================
def rotate_vector(v_rel, chi, phi):
    """
    相対速度ベクトルを散乱角(χ, φ)で回転させる
    
    Args:
        v_rel: 相対速度ベクトル (n_samples, 3) [m/s]
        chi: 散乱角 (n_samples,) [rad]
        phi: 方位角 (n_samples,) [rad]
    
    Returns:
        v_rel_new: 散乱後の相対速度ベクトル (n_samples, 3) [m/s]
    """
    n = v_rel.shape[0]
    v_rel_mag = np.linalg.norm(v_rel, axis=1, keepdims=True)
    
    # ゼロベクトルへの対処
    v_rel_mag = np.where(v_rel_mag < 1e-30, 1e-30, v_rel_mag)
    
    # 単位ベクトル
    e_z = v_rel / v_rel_mag  # (n, 3)
    
    # e_z に垂直な基底ベクトルを構築
    # 任意の垂直ベクトルを作る
    e_x = np.zeros_like(e_z)
    
    # e_z が z軸に近い場合とそうでない場合で分岐
    z_aligned = np.abs(e_z[:, 2]) > 0.9
    
    # z軸に近い場合: x軸との外積
    e_x[z_aligned, 0] = -e_z[z_aligned, 1]
    e_x[z_aligned, 1] = e_z[z_aligned, 0]
    e_x[z_aligned, 2] = 0
    
    # そうでない場合: z軸との外積
    e_x[~z_aligned, 0] = -e_z[~z_aligned, 2]
    e_x[~z_aligned, 1] = 0
    e_x[~z_aligned, 2] = e_z[~z_aligned, 0]
    
    # 正規化
    e_x_mag = np.linalg.norm(e_x, axis=1, keepdims=True)
    e_x_mag = np.where(e_x_mag < 1e-30, 1e-30, e_x_mag)
    e_x = e_x / e_x_mag
    
    # e_y = e_z × e_x
    e_y = np.cross(e_z, e_x)
    
    # 散乱後の方向
    # v' = |v| * (sin(χ)*cos(φ)*e_x + sin(χ)*sin(φ)*e_y + cos(χ)*e_z)
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


# =============================================================================
# シミュレーション実行
# =============================================================================
def run_simulation(T_neutral, T_ion, n_samples, 
                   energy_grid, sigma_elastic,
                   prob_grid, energy_grid_angle, angle_cdf):
    """
    指定された温度でモンテカルロシミュレーションを実行
    
    Args:
        T_neutral: 中性粒子温度 [eV]
        T_ion: イオン温度 [eV]
        n_samples: サンプル数
        energy_grid, sigma_elastic: 弾性散乱データ
        prob_grid, energy_grid_angle, angle_cdf: 散乱角CDFデータ
    
    Returns:
        power_coef_total: 全エネルギー輸送率係数 [m³ eV/s]
        power_coef_cx: CX成分
        power_coef_el: 弾性散乱成分
    """
    mu = M_D / 2  # 換算質量（D-D衝突）
    
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
    
    # ===== 荷電交換 (CX) の寄与 =====
    sigma_cx = get_cx_sigma(E_rel)
    # CXでは速度が完全に入れ替わる
    Delta_E_cx = E_i - E_n
    P_cx = sigma_cx * v_rel_mag * Delta_E_cx
    
    # ===== 弾性散乱 (EL) の寄与 =====
    sigma_el = get_elastic_sigma(E_rel, energy_grid, sigma_elastic)
    
    # 散乱角をサンプリング
    chi = sample_scattering_angle_vectorized(
        E_rel, prob_grid, energy_grid_angle, angle_cdf, n_samples
    )
    
    # 方位角は等方的
    phi = np.random.uniform(0, 2*np.pi, n_samples)
    
    # 散乱後の相対速度ベクトル
    v_rel_new = rotate_vector(v_rel, chi, phi)
    
    # 重心速度
    V_cm = 0.5 * (v_n + v_i)
    
    # 散乱後の中性粒子速度
    v_n_new = V_cm + 0.5 * v_rel_new
    
    # 散乱後のエネルギー [eV]
    E_n_new = 0.5 * M_D * np.sum(v_n_new**2, axis=1) / EV_TO_J
    
    Delta_E_el = E_n_new - E_n
    P_el = sigma_el * v_rel_mag * Delta_E_el
    
    # ===== 集計 =====
    power_coef_cx = np.mean(P_cx)
    power_coef_el = np.mean(P_el)
    power_coef_total = power_coef_cx + power_coef_el
    
    return power_coef_total, power_coef_cx, power_coef_el


# =============================================================================
# メイン処理
# =============================================================================
def main():
    print("=" * 60)
    print("Monte Carlo Simulation: D⁰-D⁺ Energy Transport")
    print("=" * 60)
    
    # CDFファイルを読み込む
    print("\nLoading CDF file...")
    cdf_file = "dd_00_elastic.cdf"
    try:
        energy_grid, sigma_elastic, prob_grid, energy_grid_angle, angle_cdf = \
            load_elastic_cdf(cdf_file)
        print(f"  Elastic cross-section: {len(energy_grid)} energy points")
        print(f"  Scattering angle CDF: {angle_cdf.shape} (prob × energy)")
        print(f"  Cross-section range: {sigma_elastic.min():.2e} - {sigma_elastic.max():.2e} m²")
    except Exception as e:
        print(f"Error loading CDF file: {e}")
        return
    
    # 温度スキャン
    T_neutral_arr = np.logspace(np.log10(0.1), np.log10(10.0), 20)
    
    print(f"\nRunning simulation:")
    print(f"  Ion temperature: {T_ION} eV (fixed)")
    print(f"  Neutral temperature: {T_neutral_arr[0]:.2f} - {T_neutral_arr[-1]:.2f} eV ({len(T_neutral_arr)} points)")
    print(f"  Samples per point: {N_SAMPLES}")
    
    results_total = []
    results_cx = []
    results_el = []
    
    for i, T_n in enumerate(T_neutral_arr):
        power_total, power_cx, power_el = run_simulation(
            T_n, T_ION, N_SAMPLES,
            energy_grid, sigma_elastic,
            prob_grid, energy_grid_angle, angle_cdf
        )
        results_total.append(power_total)
        results_cx.append(power_cx)
        results_el.append(power_el)
        
        print(f"  [{i+1:2d}/{len(T_neutral_arr)}] Tn = {T_n:6.3f} eV: "
              f"Total = {power_total:+.3e}, CX = {power_cx:+.3e}, EL = {power_el:+.3e}")
    
    results_total = np.array(results_total)
    results_cx = np.array(results_cx)
    results_el = np.array(results_el)
    
    # プロット
    print("\nGenerating plot...")
    
    fig, ax = plt.subplots(figsize=(10, 7))
    
    ax.plot(T_neutral_arr, results_total, 'k-', linewidth=2, label='Total')
    ax.plot(T_neutral_arr, results_cx, 'b--', linewidth=2, label='Charge Exchange (CX)')
    ax.plot(T_neutral_arr, results_el, 'r:', linewidth=2, label='Elastic Scattering (EL)')
    
    ax.axhline(y=0, color='gray', linestyle='-', linewidth=0.5)
    ax.axvline(x=T_ION, color='green', linestyle='--', linewidth=1, alpha=0.7, label=f'Ti = {T_ION} eV')
    
    ax.set_xlabel('Neutral Temperature $T_n$ [eV]', fontsize=12)
    ax.set_ylabel('Power Transfer Coefficient $\\langle \\sigma v \\Delta E \\rangle$ [m³ eV/s]', fontsize=12)
    ax.set_title('D⁰-D⁺ Energy Transport Rate Coefficient', fontsize=14)
    ax.legend(loc='best', fontsize=10)
    ax.grid(True, alpha=0.3)
    
    # 注釈を追加
    ax.annotate('Heating\n(Tn < Ti)', xy=(0.2, 0.7), xycoords='axes fraction',
                fontsize=10, ha='center', color='darkgreen')
    ax.annotate('Cooling\n(Tn > Ti)', xy=(0.8, 0.3), xycoords='axes fraction',
                fontsize=10, ha='center', color='darkred')
    
    plt.tight_layout()
    
    # 保存
    output_file = 'energy_transport_result.png'
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"Plot saved to: {output_file}")
    
    plt.show()
    
    print("\n" + "=" * 60)
    print("Simulation completed!")
    print("=" * 60)


if __name__ == "__main__":
    main()
