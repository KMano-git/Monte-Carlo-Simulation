import numpy as np
import matplotlib.pyplot as plt
import os

def sample_maxwell_velocity(temperature, mass_kg, n_samples):
    """マクスウェル分布から速度をサンプリング"""
    EV_TO_J = 1.60218e-19
    v_th = np.sqrt(temperature * EV_TO_J / mass_kg)
    velocity = np.random.normal(0, v_th, size=(n_samples, 3))
    return velocity

def main():
    # 物理定数
    AMU = 1.66054e-27
    M_D = 2.014 * AMU
    EV_TO_J = 1.60218e-19
    
    # シミュレーション設定
    T_neutral = 2.0  # 中性粒子温度 [eV]
    T_ion = 2.0      # イオン温度 [eV]
    N_SAMPLES = int(2e6)
    
    print("============================================================")
    print("Charge Exchange (CX) Scattering Angle Validation")
    print("============================================================")
    print("  Simulation implementation logic:")
    print("    collision_cx: p%vx = vx_i")
    print("  This script replicates this logic to compute scattering angles.")
    
    # 1. 衝突前の中性粒子速度 (v_n) と イオン速度 (v_i) をサンプリング
    v_n = sample_maxwell_velocity(T_neutral, M_D, N_SAMPLES)
    # 背景プラズマの流速が0と仮定したマクスウェル分布
    v_i = sample_maxwell_velocity(T_ion, M_D, N_SAMPLES)
    
    # CXの判定や重み（断面積×相対速度）は今回は角度分布自体のデモなので
    # 十分大きなサンプル数の一様イベントとして扱います。
    # 厳密には反応率σvで重み付けされますが、T_n=T_i=2.0eV で流速0であれば、
    # 等方性（Lab系での無相関性）は維持されます。
    
    # 相対速度
    v_rel = v_n - v_i
    v_rel_mag = np.linalg.norm(v_rel, axis=1)
    
    # 2. CX処理（Fortranと全く同じロジック）
    # 新しい中性粒子は、衝突したイオンの速度をそのまま受け継ぐ
    v_n_new = v_i
    
    # ---------------------------------------------------------
    # 室系 (Lab Frame) での散乱角の計算
    # ---------------------------------------------------------
    v_n_mag = np.linalg.norm(v_n, axis=1)
    v_n_new_mag = np.linalg.norm(v_n_new, axis=1)
    
    # ゼロ割りを防ぐ
    valid = (v_n_mag > 0) & (v_n_new_mag > 0)
    
    cos_chi_lab = np.sum(v_n[valid] * v_n_new[valid], axis=1) / (v_n_mag[valid] * v_n_new_mag[valid])
    cos_chi_lab = np.clip(cos_chi_lab, -1.0, 1.0)
    chi_lab_deg = np.degrees(np.arccos(cos_chi_lab))
    
    # ---------------------------------------------------------
    # 重心系 (CM Frame) での散乱角の計算
    # ---------------------------------------------------------
    # 重心速度 V_cm = (m_n*v_n + m_i*v_i) / (m_n + m_i)
    # 質量は等しい(M_D)ので単なる平均
    V_cm = 0.5 * (v_n + v_i)
    
    # 重心系での衝突前の速度
    v_n_cm = v_n - V_cm       # = 0.5 * (v_n - v_i)
    
    # 重心系での衝突後の速度 (v_n_new = v_i なので)
    v_n_new_cm = v_n_new - V_cm # = v_i - 0.5 * (v_n + v_i) = -0.5 * (v_n - v_i) = - v_n_cm
    
    v_n_cm_mag = np.linalg.norm(v_n_cm, axis=1)
    v_n_new_cm_mag = np.linalg.norm(v_n_new_cm, axis=1)
    
    valid_cm = (v_n_cm_mag > 0) & (v_n_new_cm_mag > 0)
    
    cos_chi_cm = np.sum(v_n_cm[valid_cm] * v_n_new_cm[valid_cm], axis=1) / (v_n_cm_mag[valid_cm] * v_n_new_cm_mag[valid_cm])
    cos_chi_cm = np.clip(cos_chi_cm, -1.0, 1.0)
    chi_cm_deg = np.degrees(np.arccos(cos_chi_cm))

    # 3. プロットの作成
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # ===== Lab Frame Plot =====
    axes[0].hist(
        chi_lab_deg, bins=90, range=(0, 180), density=True,
        color='blue', alpha=0.6, edgecolor='darkblue', label='Simulation (CX Lab Frame)'
    )
    
    # 理論的な等方散乱の確率密度関数: P(θ) = (1/2) * sin(θ)  [rad基準] -> 度数法に変換
    x_angles = np.linspace(0, 180, 200)
    # P_deg(θ) d(θ_deg) = 1/2 sin(θ) d(θ_rad) = 1/2 sin(θ) * (π/180) d(θ_deg)
    iso_pdf = 0.5 * np.sin(np.radians(x_angles)) * (np.pi / 180.0)
    axes[0].plot(x_angles, iso_pdf, 'k--', linewidth=2, label='Isotropic Theory (sinθ)')
    
    axes[0].set_title('Scattering Angle in LAB Frame')
    axes[0].set_xlabel('Scattering Angle [degrees]')
    axes[0].set_ylabel('Probability Density')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    # ===== CM Frame Plot =====
    # 理論的には全てがちょうど180度になるはずですが、ヒストグラムで確認
    axes[1].hist(
        chi_cm_deg, bins=90, range=(0, 180), density=True,
        color='red', alpha=0.6, edgecolor='darkred', label='Simulation (CX CM Frame)'
    )
    
    axes[1].axvline(180.0, color='k', linestyle='--', linewidth=2, label='Expected: 180°')
    
    axes[1].set_title('Scattering Angle in CM Frame')
    axes[1].set_xlabel('Scattering Angle [degrees]')
    axes[1].set_ylabel('Probability Density')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    output_filename = 'cx_scattering_angle_validation.png'
    plt.savefig(output_filename, dpi=200)
    print(f"\nPlot saved to: {output_filename}")
    
    # 統計出力
    print("\n===== Statistical Results =====")
    print(f"LAB Frame Mean Angle: {np.mean(chi_lab_deg):.2f} degrees (Isotropic expected ~ 90.0)")
    print(f"CM Frame Mean Angle:  {np.mean(chi_cm_deg):.2f} degrees (Backscatter expected = 180.0)")

if __name__ == '__main__':
    main()
