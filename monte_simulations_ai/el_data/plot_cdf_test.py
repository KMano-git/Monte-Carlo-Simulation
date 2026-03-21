"""
dd_00_elastic_pure_el_angle_fixed.cdf の各変数を図示するスクリプト

CDFファイル構造 (xs_data_tab 内の格納順):
  [0] cross_section    : rank=1, 101点,       E [eV/amu]
  [1] reaction_rate    : rank=2, 51x51,       E [eV/amu] x T [eV/amu]
  [2] I_1_0            : rank=2, 51x51,       E [eV/amu] x T [eV/amu]
  [3] I_1_1*up         : rank=2, 51x51,       E [eV/amu] x T [eV/amu]
  [4] I_1_2*up^2       : rank=2, 51x51,       E [eV/amu] x T [eV/amu]
  [5] scattering_angle : rank=2, 251x51,      random_num x E [eV/amu]
  [6] sigv_max         : rank=0, scalar
  [7] angle_min        : rank=0, scalar
"""
import re
import numpy as np
import matplotlib.pyplot as plt
import os

# =====================================================================
# CDFファイルのパース
# =====================================================================
CDF_PATH = "dd_00_elastic_pure_el_angle_fixed.cdf"
OUT_DIR  = os.path.dirname(os.path.abspath(__file__))

def parse_array(text, name):
    """CDL テキストからカンマ区切りの数値配列を抽出する"""
    pat = re.compile(rf'\b{re.escape(name)}\s*=\s*(.*?)\s*;', re.DOTALL)
    m = pat.search(text)
    if not m:
        raise ValueError(f"'{name}' not found in CDF")
    raw = m.group(1).replace('\n', ' ')
    vals = [x.strip() for x in raw.split(',') if x.strip()]
    return vals

def parse_int_array(text, name):
    return [int(v) for v in parse_array(text, name)]

def parse_float_array(text, name):
    return [float(v) for v in parse_array(text, name)]

with open(CDF_PATH, 'r') as f:
    cdf_text = f.read()

xs_rank      = parse_int_array(cdf_text, 'xs_rank')
xs_data_base = parse_int_array(cdf_text, 'xs_data_base')
xs_data_inc  = parse_int_array(cdf_text, 'xs_data_inc')
xs_data_tab  = parse_float_array(cdf_text, 'xs_data_tab')

# xs_tab_index は (20, 3) の2次元配列
tab_flat = parse_int_array(cdf_text, 'xs_tab_index')
xs_tab_index = [tab_flat[i*3:(i+1)*3] for i in range(20)]

# xs_min, xs_max は (20, 3) の2次元配列
min_flat = parse_float_array(cdf_text, 'xs_min')
max_flat = parse_float_array(cdf_text, 'xs_max')
xs_min = [min_flat[i*3:(i+1)*3] for i in range(20)]
xs_max = [max_flat[i*3:(i+1)*3] for i in range(20)]

# xs_spacing は (20, 4) の2次元配列（文字列）
spacing_raw = parse_array(cdf_text, 'xs_spacing')
xs_spacing = [spacing_raw[i*4:(i+1)*4] for i in range(20)]

# 変数メタデータ
var_info = [
    {"name": "cross_section",    "unit": "cm^2",    "idx": 0},
    {"name": "reaction_rate",    "unit": "cm^3/s",   "idx": 1},
    {"name": "I_1_0",            "unit": "cm^3/s",   "idx": 2},
    {"name": "I_1_1*up",         "unit": "cm^4/s^2", "idx": 3},
    {"name": "I_1_2*up^2",       "unit": "cm^5/s^3", "idx": 4},
    {"name": "scattering_angle", "unit": "radians",  "idx": 5},
    {"name": "sigv_max",         "unit": "cm^3/s",   "idx": 6},
    {"name": "angle_min",        "unit": "radians",  "idx": 7},
]

def get_data(idx):
    """xs_data_tab から idx 番目の変数のデータを取り出す"""
    base = xs_data_base[idx]
    size = xs_data_inc[idx]
    return np.array(xs_data_tab[base:base+size])

def make_axis(n, vmin, vmax, spacing):
    """グリッド軸を生成する (log or linear)"""
    spacing = spacing.replace('"', '').strip()
    if spacing == 'log' and vmin > 0 and vmax > 0:
        return np.logspace(np.log10(vmin), np.log10(vmax), n)
    else:
        return np.linspace(vmin, vmax, n)

# =====================================================================
# 1. cross_section (1D) : σ(E)
# =====================================================================
print("[1/6] cross_section ...")
idx = 0
data = get_data(idx)
n1 = xs_tab_index[idx][0]  # 101
E_cs = make_axis(n1, xs_min[idx][0], xs_max[idx][0], xs_spacing[idx][0])

fig, ax = plt.subplots(figsize=(8, 5))
ax.loglog(E_cs, data, 'b-', lw=1.5)
ax.set_ylim(1e-15, 1e-12)
ax.set_xlabel("E [eV/amu]")
ax.set_ylabel("cross_section [cm²]")
ax.set_title("cross_section  σ(E)")
ax.grid(True, which='both', alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "cross_section.png"), dpi=150)
plt.close(fig)

# =====================================================================
# 2D → 代表温度での線プロット用の共通関数
# =====================================================================
def plot_2d_lines(idx, title, fname, ylabel):
    """2D データを代表温度ごとの線プロットとして描画する"""
    data = get_data(idx)
    n1 = xs_tab_index[idx][0]
    n2 = xs_tab_index[idx][1]
    Z = data.reshape(n2, n1)  # shape = (nT, nE)

    E = make_axis(n1, xs_min[idx][0], xs_max[idx][0], xs_spacing[idx][0])
    T = make_axis(n2, xs_min[idx][1], xs_max[idx][1], xs_spacing[idx][1])

    # 5点を等間隔インデックスで選ぶ（被りなし）
    t_indices = np.linspace(0, n2 - 1, 5, dtype=int)

    fig, ax = plt.subplots(figsize=(9, 6))
    for j in t_indices:
        T_actual = T[j]
        ax.loglog(E, np.abs(Z[j, :]), '-', lw=1.5,
                  label=f"T = {T_actual:.4g} eV/amu")

    ax.set_xlabel("E [eV/amu]")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(fontsize=9)
    ax.grid(True, which='both', alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, fname), dpi=150)
    plt.close(fig)

# =====================================================================
# 2. reaction_rate
# =====================================================================
print("[2/6] reaction_rate ...")
plot_2d_lines(1, r"reaction_rate  $\langle\sigma v\rangle$(E)  at various T",
              "reaction_rate.png", "reaction_rate [cm³/s]")

# =====================================================================
# 3. I_1_0
# =====================================================================
print("[3/6] I_1_0 ...")
plot_2d_lines(2, r"$I_{1,0}$(E)  at various T",
              "I_1_0.png", r"$I_{1,0}$ [cm³/s]")

# =====================================================================
# 4. I_1_1*up
# =====================================================================
print("[4/6] I_1_1*up ...")
plot_2d_lines(3, r"$I_{1,1} \cdot u_\parallel$(E)  at various T",
              "I_1_1_up.png", r"$I_{1,1} \cdot u_\parallel$ [cm⁴/s²]")

# =====================================================================
# 5. I_1_2*up^2
# =====================================================================
print("[5/6] I_1_2*up^2 ...")
plot_2d_lines(4, r"$I_{1,2} \cdot u_\parallel^2$(E)  at various T",
              "I_1_2_up2.png", r"$I_{1,2} \cdot u_\parallel^2$ [cm⁵/s³]")

# =====================================================================
# 6. scattering_angle — 代表エネルギーでの線プロット
# =====================================================================
print("[6/6] scattering_angle ...")
idx = 5
data = get_data(idx)
n1 = xs_tab_index[idx][0]  # 251
n2 = xs_tab_index[idx][1]  # 51
Z = data.reshape(n2, n1)   # shape = (nE, n_rand)

x_rand = make_axis(n1, xs_min[idx][0], xs_max[idx][0], xs_spacing[idx][0])
E_scat = make_axis(n2, xs_min[idx][1], xs_max[idx][1], xs_spacing[idx][1])

# 5点を等間隔インデックスで選ぶ（被りなし）
e_indices = np.linspace(0, n2 - 1, 5, dtype=int)

fig, ax = plt.subplots(figsize=(9, 6))
for j in e_indices:
    E_actual = E_scat[j]
    ax.plot(x_rand, Z[j, :], '-', lw=1.5,
            label=f"E = {E_actual:.4g} eV/amu")

ax.set_xlabel("1st random number (CDF)")
ax.set_ylabel("scattering angle θ [rad]")
ax.set_ylim(0, np.pi)
ax.set_title("scattering_angle  θ(rand)  at various E")
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "scattering_angle.png"), dpi=150)
plt.close(fig)

# =====================================================================
# スカラ値の出力
# =====================================================================
sigv_max = xs_data_tab[xs_data_base[6]]
angle_min = xs_data_tab[xs_data_base[7]]
print(f"\nsigv_max  = {sigv_max:.6e} cm^3/s")
print(f"angle_min = {angle_min:.6e} rad")

print(f"\n全図を {OUT_DIR} に保存しました。")
