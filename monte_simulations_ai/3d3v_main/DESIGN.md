# 3D3V Monte Carlo Simulation ver 1.0.0 設計書

## 1. 目的

`3d3v_ver1.0.0` は、`../3d3v_lookup_check` を基準に、外部発表に使える 3D3V 非アナログ Monte Carlo シミュレーションとして整理する。

当面の主目的は、弾性散乱断面積データの違いが、イオンエネルギー、運動量、粒子ソースなどに与える影響を比較できるようにすることである。

本設計では、既存コードから変更する主対象を EL 項に限定する。EI と CX の基本アルゴリズムは、現行 `3d3v_lookup_check` および `src_read/monte` の考え方を維持する。

## 2. 基準コード

実装の基準は次とする。

- 基準ディレクトリ: `../3d3v_lookup_check`
- 参照する元コード:
  - `../../src_read/monte/ntscor_cl.f`
  - `../../src_read/monte/ntscor_tr.f`
  - `../../src_read/monte/ntscrg.f`
  - `../../src_read/monte/ntsumsr.f`
  - `../../src_read/elcol/ntvel_el.f`
  - `../../src_read/elcol/ntsctr_el.f`

`3d3v_ver1.0.0` の既存ファイル状態は前提にしない。README の方針と `3d3v_lookup_check` の動作をもとに新しく形にする。

## 3. 方針

### 3.1 Estimator 名称

`3d3v_lookup_check` での旧名称と、新コードでの名称を次のように対応させる。

```text
旧 A(EL)        -> 新 CL
旧 CL           -> 新 CL(avg)
旧 TR           -> 新 TR
旧 CLInner      -> 削除
旧 TRInner      -> 削除
旧 TRPretab     -> compare_estimator = .true. の場合のみ比較用に出力
```

最終的な通常出力は、次の 3 estimator を主系列とする。

```text
CL
CL(avg)
TR
```

`compare_estimator = .true.` の場合だけ、`TR(pretab)` などの比較用列を追加する。

### 3.2 変更対象

今回の物理モデル変更対象は EL 項とする。

EI と CX は、CL、CL(avg)、TR の全 estimator で同じ定義を使う。CL と CL(avg) の差は EL 項の score 方法だけである。

### 3.3 削除する機能

次は ver1.0.0 では削除する。

- inner Monte Carlo estimator
- `tl_el_inner_samples`
- `use_isotropic`
- 等方散乱分岐

EL 散乱角は常に CDF データからサンプリングする。

## 4. 物理量と符号規約

`src_read/monte` では plasma 側の source を次の 4 種類で扱っている。

```text
Sn : plasma ion particle source
Sp : plasma ion momentum source
We : electron energy source
Wi : ion energy source
```

新 3D3V コードでもこの分類を採用する。

符号は plasma 側を基準にする。

```text
正の Sn : plasma ion が増える
正の Sp : plasma ion が運動量を受け取る
正の We : electron fluid がエネルギーを受け取る
正の Wi : ion fluid がエネルギーを受け取る
```

ただし `src_read/monte` の EI や解離反応では、電子エネルギー損失を `We < 0` として記録している。したがって EI では `We = - electron_loss` になる。

本コードは 3D3V なので、`Sp` は平行成分だけでなく 3 成分を持つ。

```text
Sp = (Spx, Spy, Spz)
```

磁力線方向成分が必要な場合は、後処理で

```text
Sp_parallel = Sp dot b
```

として求める。

## 5. 共通記号

```text
w        : Monte Carlo particle weight
n_i      : background ion density [m^-3]
n_e      : electron density [m^-3]
T_i      : ion temperature [eV]
T_e      : electron temperature [eV]
m_D      : deuterium mass [kg]
e        : elementary charge [J/eV]
v0       : neutral velocity before reaction [m/s]
v0'      : neutral velocity after reaction [m/s]
vi       : sampled background ion velocity [m/s]
u_i      : background ion flow velocity [m/s]
E0       : 0.5 * m_D * |v0|^2 [J]
E0'      : 0.5 * m_D * |v0'|^2 [J]
Ei_mean  : 1.5 * T_i * e + 0.5 * m_D * |u_i|^2 [J]
R_a      : electron impact ionization rate [s^-1]
R_cx     : charge-exchange rate [s^-1]
R_el     : elastic collision rate [s^-1]
R_s      : R_cx + R_el [s^-1]
tau_eff  : integral w(t) dt over a flight segment [s]
```

`tau_eff` は非アナログ EI による weight 減衰を含む有効飛行時間である。現行 `3d3v_lookup_check` と同じく、

```text
tau_eff = w_old * (1 - exp(-R_a * dt)) / R_a    if R_a > 0
tau_eff = w_old * dt                            if R_a = 0
```

を使う。

## 6. 反応別 source 定義

### 6.1 EI: e + D0 -> D+ + 2e

EI は非アナログ法で neutral weight を連続的に減衰させる。

`src_read/monte/ntscor_cl.f` では、EI の粒子数を

```text
dN = weitb - weit
```

として扱う。新コードでは estimator ごとに次を使う。

CL と CL(avg) では、accepted scatter collision ごとに

```text
dN_ei = w * (R_a / R_s)
```

を score する。これは accepted scatter collision の発生頻度が `R_s` であるため、時間平均では `w * R_a` に一致する。

TR では flight segment ごとに

```text
dN_ei = R_a * tau_eff
```

を score する。

EI の source は次とする。

```text
Sn_i += dN_ei

Sp_i += dN_ei * m_D * v0

We   += -dN_ei * E_e_loss

Wi   +=  dN_ei * E0
```

ここで `E_e_loss` は electron impact ionization による電子エネルギー損失である。現行 3D3V の簡略モデルで詳細な electron loss table を使わない段階では、少なくとも電離ポテンシャルを使う。

```text
E_e_loss = E_IONIZE_THRESHOLD * e
```

将来 `src_read/monte` の `elrc` に相当する電子エネルギー損失データを導入した場合は、

```text
E_e_loss = elrc(E, T_e, n_e) * e
```

に置き換える。

既存のイオンエネルギーだけの score と対応させる場合、`Wi` の EI 項は中性粒子が ion fluid に持ち込む運動エネルギーである。

```text
Wi_EI = dN_ei * 0.5 * m_D * |v0|^2
```

### 6.2 CX: D0 + D+ -> D+ + D0

CX は、飛んでいる neutral が plasma ion になり、背景 ion が新しい neutral になる交換反応である。

CL と CL(avg) では、accepted scatter collision ごとに

```text
dN_cx = w * (R_cx / R_s)
```

を使う。

TR では、背景 ion 速度を sampling する実装では

```text
dN_cx = tau_eff * n_i * sigma_cx(E_rel) * v_rel
```

を使う。

ここで

```text
v_rel = |v0 - vi|
E_rel = 0.25 * m_D * v_rel^2 / e
```

である。

CX の particle source は、単一 D species では net 0 である。ただし保存則確認のため、生成と消滅を別々に保持する。

```text
Sn_i_plus  += dN_cx      old neutral -> ion
Sn_i_minus += dN_cx      background ion -> neutral
Sn_i_net   += 0
```

運動量 source は、3D3V では 3 成分で保持する。

sampled ion velocity を使う場合は、

```text
Sp_i += dN_cx * m_D * v0
Sp_i -= dN_cx * m_D * vi
```

`src_read/monte` の fluid source と同じ平均量を使う場合は、背景 ion 速度を flow 平均に置き換える。

```text
Sp_i += dN_cx * m_D * v0
Sp_i -= dN_cx * m_D * u_i
```

イオンエネルギー source も同様に、sampling 版と fluid 平均版を区別する。

sampling 版:

```text
Wi += dN_cx * (0.5 * m_D * |v0|^2)
Wi -= dN_cx * (0.5 * m_D * |vi|^2)
```

fluid 平均版:

```text
Wi += dN_cx * E0
Wi -= dN_cx * Ei_mean
```

`src_read/monte/ntscor_cl.f` と `ntscor_tr.f` では、CX の ion energy source は次の形で書かれている。

```text
Wi_CX = dN_cx * E0
      - dN_cx * (1.5 * T_i * e + 0.5 * m_D * |u_i|^2)
```

本 3D3V_ver1.0.0 の通常 source 出力では、EI/CX を現行 `3d3v_lookup_check` から変更しない方針を優先し、sampled ion velocity 版を標準とする。

したがって、CL と CL(avg) は accepted collision 判定で実際に sample した `vi` を使い、TR は CX rate 評価のために sample した `vi` を同じ source score に使う。

```text
Wi_CX_sample = dN_cx * (E0 - 0.5 * m_D * |vi|^2)
```

fluid 平均版は `src_read/monte` との比較用の別 diagnostic として追加できるが、EL だけを変更する検証では標準出力に混ぜない。

### 6.3 EL: D0 + D+ -> D0 + D+

EL では粒子種は変わらない。

```text
Sn = 0
We = 0
```

neutral の運動量・エネルギー変化を

```text
Delta_p0 = m_D * (v0' - v0)
Delta_E0 = E0' - E0
```

とする。

plasma ion 側の source は反作用なので、

```text
Sp_i += -w * Delta_p0
Wi   += -w * Delta_E0
```

である。

#### CL の EL

新 `CL` は旧 `A(EL)` に対応する。実際に EL collision が選ばれた場合だけ、実現値を score する。

```text
if actual channel is EL:
    Delta_p0 = m_D * (v0_after - v0_before)
    Delta_E0 = 0.5 * m_D * |v0_after|^2
             - 0.5 * m_D * |v0_before|^2

    Sp_EL_CL += -w * Delta_p0
    Wi_EL_CL += -w * Delta_E0
```

この estimator は EL の実衝突サンプリング揺らぎをそのまま持つ。

#### CL(avg) の EL

新 `CL(avg)` は旧 `CL` に対応する。accepted scatter collision ごとに、EL の期待寄与を score する。

```text
dN_el = w * (R_el / R_s)

Sp_EL_CLavg += -dN_el * <Delta_p0>_EL
Wi_EL_CLavg += -dN_el * <Delta_E0>_EL
```

`<Delta_p0>_EL` と `<Delta_E0>_EL` は、CDF の EL lookup kernel から求める。現行 `3d3v_lookup_check` の energy 実装では、`compute_el_lookup_rate_and_score` により

```text
el_rate            = n_i * < -Delta_E0 * sigma_el * v_rel >
el_reaction_rate   = n_i * < sigma_el * v_rel >
el_collision_score = el_rate / el_reaction_rate
```

を作っている。`CL(avg)` ではこの `el_collision_score` を EL 1 collision あたりの期待 ion energy source として使う。

運動量も同じ考え方で、EL lookup kernel に

```text
el_momentum_rate[3] = n_i * < -Delta_p0[3] * sigma_el * v_rel >
el_momentum_score[3] = el_momentum_rate[3] / el_reaction_rate
```

を追加する。

#### TR の EL

TR は flight segment ごとに、単位時間 source を `tau_eff` に掛ける。

```text
Sp_EL_TR += el_momentum_rate[3] * tau_eff
Wi_EL_TR += el_energy_rate      * tau_eff
```

`src_read/elcol/ntsctr_el.f` では、EL の track-length estimator は次の量を計算している。

```text
dtr_mvp = <Delta(m v_parallel) * sigma_el * v_rel>
dtr_eng = <Delta(E0)           * sigma_el * v_rel>
```

source 側では反作用として

```text
Sp_parallel += -N0 * n_i * dtr_mvp
Wi          += -N0 * n_i * dtr_eng
```

を使う。

新 3D3V では、平行成分だけでなく 3 成分を保持する。

```text
dtr_mv[3] = <Delta_p0[3] * sigma_el * v_rel>

Sp_EL_TR[3] = -n_i * dtr_mv[3] * tau_eff
Wi_EL_TR    = -n_i * dtr_eng  * tau_eff
```

現行 CDF の `I_1_0`, `I_1_1*up`, `I_1_2*up^2` から energy rate を作る式は `3d3v_lookup_check` の `compute_el_lookup_rate_and_score` を維持する。運動量 rate については `src_read/elcol/ntsctr_el.f` の

```text
sp0 = I_1_1 / |ut| - 0.5 * va2 / |ut|^2 * I_1_0
dtr_mv = -ut * reduced_mass * sp0
```

を 3 成分に拡張して使う。

ここで

```text
ut  = v0 - u_i
va2 = 2 * T_i * e / m_D
```

である。

## 7. Estimator ごとの実装アルゴリズム

### 7.1 CL

`CL` は EI/CX を旧 CL のまま、EL を actual collision score にする。

accepted scatter collision ごとに、衝突前状態 `p_before` を使って EI/CX を score する。

```text
R_cx = n_i * sigma_cx(E_rel) * v_rel
R_el = n_i * sigma_el(E_rel) * v_rel
R_s  = R_cx + R_el
R_a  = n_e * <sigma v>_ion(T_e)

score_ei_common(dN = w * R_a  / R_s)
score_cx_common(dN = w * R_cx / R_s)
```

その後、actual channel が EL の場合だけ、

```text
score_el_actual(p_before, p_after)
```

を呼ぶ。

actual channel が CX の場合でも EI/CX common score は行う。これは旧 CL の扱いと同じである。

### 7.2 CL(avg)

`CL(avg)` も EI/CX は `CL` と同じ common score を使う。

accepted scatter collision ごとに、

```text
score_ei_common(dN = w * R_a  / R_s)
score_cx_common(dN = w * R_cx / R_s)
score_el_average(dN = w * R_el / R_s)
```

を呼ぶ。

`score_el_average` は CDF lookup から EL 1 collision あたりの期待寄与を得る。

### 7.3 TR

TR は flight segment の flush 時に score する。

```text
tau_eff = integral w(t) dt
```

EI:

```text
dN_ei = R_a * tau_eff
score_ei_common(dN_ei)
```

CX:

```text
dN_cx_rate = n_i * <sigma_cx * v_rel>
dN_cx      = dN_cx_rate * tau_eff
score_cx_common(dN_cx)
```

現行実装では、`3d3v_lookup_check` と同じく背景 ion velocity を sample して CX TR を評価し、その同じ sampled ion velocity で source も score する。

EL:

```text
score_el_track_length(tau_eff)
```

を呼び、CDF lookup から `Sp` と `Wi` の rate を得る。

## 8. データ構造

### 8.1 反応種

ver1.0.0 の最初の対象は D0 の次の 3 反応とする。

```text
EI
CX
EL
TOTAL
```

### 8.2 Source データ

source は estimator ごと、反応ごとに保持する。

```fortran
type :: source_terms_t
   real(dp) :: sn_plus
   real(dp) :: sn_minus
   real(dp) :: sn_net
   real(dp) :: sp(3)
   real(dp) :: we
   real(dp) :: wi
end type source_terms_t
```

`sn_net` は

```text
sn_net = sn_plus - sn_minus
```

である。

estimator 全体は次のように持つ。

```fortran
type :: estimator_score_t
   type(source_terms_t) :: ei
   type(source_terms_t) :: cx
   type(source_terms_t) :: el
   type(source_terms_t) :: total
end type estimator_score_t

type :: score_data
   type(estimator_score_t) :: cl
   type(estimator_score_t) :: cl_avg
   type(estimator_score_t) :: tr
   type(estimator_score_t) :: tr_pretab
end type score_data
```

`tr_pretab` は `compare_estimator = .true.` のときだけ更新し、通常出力では空でよい。

### 8.3 統計データ

標準偏差と信頼区間のため、step score または block score を統計単位として蓄積する。

```fortran
type :: running_stats_t
   integer :: n
   real(dp) :: mean
   real(dp) :: m2
end type running_stats_t
```

Welford 法で平均と分散を更新する。

```text
variance = m2 / (n - 1)
stderr   = sqrt(variance / n)
ci95     = 1.96 * stderr
```

標準の統計単位は timestep とする。粒子ごとの統計が必要な場合は、粒子 local score を同じ構造に流し込めるようにする。

## 9. 入力パラメータ

`run/input.nml` の `&simulation` は次を基本にする。

```fortran
&simulation
    n_particles = 10000
    n_steps     = 100
    dt          = 1.0d-9
    seed        = 12345

    enable_cx = .true.
    enable_el = .true.
    enable_ei = .true.

    compare_estimator = .false.

    weight_min = 1.0d-10

    cdf_file         = 'dd_00_elastic.cdf'
    tl_el_table_file = 'tl_el_table.dat'

    output_sources   = 'sources.csv'
    output_stats     = 'source_stats.csv'
    output_balance   = 'balance.csv'
    output_hist      = 'energy_hist.csv'
/
```

削除する入力:

```text
tl_el_inner_samples
use_isotropic
```

## 10. 出力仕様

### 10.1 `sources.csv`

各 timestep の source を出力する。

最小列は次の通り。

```text
time[s]
estimator
reaction
Sn_plus[m^-3 s^-1]
Sn_minus[m^-3 s^-1]
Sn_net[m^-3 s^-1]
Spx[N m^-3]
Spy[N m^-3]
Spz[N m^-3]
We[W m^-3]
Wi[W m^-3]
```

内部 score から物理単位への変換は、現行 `ntscrg.csv` と同じ考え方で行う。

```text
source_density = score * n_init / (N_particles * dt)
```

quantity ごとの単位は次である。

```text
Sn score : dimensionless weight
Sp score : kg m/s
We score : J
Wi score : J
```

したがって出力単位は

```text
Sn : m^-3 s^-1
Sp : kg m^-2 s^-2 = N m^-3
We : W m^-3
Wi : W m^-3
```

### 10.2 `source_stats.csv`

step score または block score から統計値を出す。

```text
estimator
reaction
quantity
mean
stddev
stderr
ci95_low
ci95_high
n_blocks
```

### 10.3 `balance.csv`

保存則確認用の収支を出す。

```text
quantity
initial
final
plasma_source_integral
neutral_change
residual
relative_residual
```

対象は次。

```text
particle
momentum_x
momentum_y
momentum_z
energy
```

CX の粒子 source は net 0 になるため、`Sn_plus` と `Sn_minus` の両方を確認する。

### 10.4 既存互換出力

既存比較用に、簡易 summary として次を残してよい。

```text
CL_Wi_ei, CL_Wi_cx, CL_Wi_el, CL_Wi_total
CLavg_Wi_ei, CLavg_Wi_cx, CLavg_Wi_el, CLavg_Wi_total
TR_Wi_ei, TR_Wi_cx, TR_Wi_el, TR_Wi_total
```

ただし主出力は `sources.csv` とする。

## 11. モジュール構成

`3d3v_lookup_check` から移植したあと、次の構成へ整理する。

```text
code/constants.f90
code/data_types.f90
code/random_utils.f90
code/cdf_reader.f90
code/cross_sections.f90
code/dynamics.f90
code/scoring.f90
code/source_terms.f90
code/statistics.f90
code/io.f90
code/main.f90
```

役割:

```text
data_types.f90
  particle_t, sim_params, plasma_params, score_data を定義する。

main.f90
  timestep と particle loop を管理する。
  OpenMP 有効時は timestep 内の particle loop だけを並列化する。
  score_data と deltaE histogram は thread-local に蓄積し、step 終了後に合算する。

dynamics.f90
  粒子移動、CX collision、EL collision、weight 更新を扱う。
  use_isotropic 分岐は持たない。

source_terms.f90
  EI/CX/EL の Sn/Sp/We/Wi 式を実装する。
  EI/CX common source はここに集約する。

scoring.f90
  CL, CL(avg), TR の estimator 呼び出しを管理する。
  物理式の細部は source_terms.f90 に寄せる。

statistics.f90
  mean, variance, stderr, confidence interval を扱う。

io.f90
  input.nml 読み込み、CSV 出力、summary 出力を扱う。
```

## 12. 実装フェーズ

### Phase 1: 基準コードの移植

`3d3v_lookup_check` から必要な source と run 入力を移植する。

持ち込む:

```text
Makefile
code/*.f90
run/input.nml
run/dd_00_elastic.cdf
```

持ち込まない:

```text
run/monte_carlo_3d3v_natl
run/build_tl_el_table
run/*.png
run/debug_output.txt
run/__pycache__
```

この段階ではまず build と run を通す。

### Phase 2: 不要機能の削除

次を削除する。

```text
score_inner_multi_collision_el
score_inner_multi_track_length
tl_el_sample_avg
tl_el_sample_stats
tl_el_inner_samples
use_isotropic
```

EL collision は常に CDF scattering angle を使う。

### Phase 3: Estimator 名称整理

score 構造を次へ置き換える。

```text
a_el       -> cl.el
cl_*       -> cl_avg.*
tr_*       -> tr.*
tr_pretab  -> tr_pretab.*
```

ただし EI/CX は `cl` と `cl_avg` の両方で同じ common source を使う。

### Phase 4: Source term 実装

`source_terms.f90` を追加し、次を実装する。

```text
add_ei_source(...)
add_cx_source(...)
add_el_actual_source(...)
add_el_average_source(...)
add_el_track_length_source(...)
```

EI/CX は CL と CL(avg) で同じ関数を呼ぶ。

### Phase 5: EL lookup の拡張

現行 energy lookup に加えて、運動量 3 成分の rate を出す。

```text
compute_el_lookup_sources(...)
  -> el_energy_rate
  -> el_momentum_rate(3)
  -> el_reaction_rate
  -> el_energy_per_collision
  -> el_momentum_per_collision(3)
```

energy rate は現行式を維持する。momentum rate は `src_read/elcol/ntsctr_el.f` の `dtr_mvx`, `dtr_mvy`, `dtr_mvz` に相当する式を実装する。

### Phase 6: 出力と統計

`sources.csv`, `source_stats.csv`, `balance.csv` を追加する。

既存の `ntscrg.csv` 形式は、互換 summary として残すか、`sources.csv` から後処理で作る。

### Phase 7: 検証

検証項目:

```text
1. compare_estimator = .false. で CL, CL(avg), TR のみ出る。
2. compare_estimator = .true. で TR(pretab) が追加される。
3. EI/CX の既存 Wi 出力が 3d3v_lookup_check と一致する。
4. EL の CL は actual EL collision の揺らぎを持つ。
5. EL の CL(avg) と TR は平均的に近い。
6. CX の Sn_net が単一 species で 0 に近い。
7. 保存則 residual が十分小さい。
8. CDF ファイルを変えたとき、EL source の差が観察できる。
```

## 13. 注意点

### 13.1 EI の電子エネルギー

現在の 3D3V コードでは、EI の電子エネルギー損失は詳細に扱っていない。`src_read/monte` では `elrc` が electron energy loss を表す。

ver1.0.0 の最初の実装では、

```text
We_EI = -dN_ei * E_IONIZE_THRESHOLD * e
```

を標準にする。より正確な電子エネルギー損失モデルは別フェーズで追加する。

### 13.2 CX の energy source

現行 `3d3v_lookup_check` の CX energy score は sampled ion velocity を使う。一方、`src_read/monte` の source 出力は Maxwell 平均の ion thermal energy と flow energyを使う。

本実装では、旧CLとの連続性を優先して sampled 版を標準にする。`src_read/monte` 準拠の fluid 平均版は、必要になった時点で optional diagnostic として追加する。

### 13.3 EL の運動量

`src_read/monte` は平行運動量 `Sp_parallel` を主に出力する。新 3D3V では 3 成分 `Spx, Spy, Spz` を主出力にする。

平行成分が必要な場合は、背景磁場方向 `b` を入力に追加し、後処理または出力時に

```text
Sp_parallel = Spx * bx + Spy * by + Spz * bz
```

を計算する。

### 13.4 compare_estimator

`TR(pretab)` は mainline ではない。通常の物理結果としては `CL`, `CL(avg)`, `TR` を見る。`TR(pretab)` は CDF lookup TR との実装比較に限って使う。

## 14. 完了条件

ver1.0.0 の初期完成条件は次とする。

```text
1. make で build できる。
2. run/input.nml から実行できる。
3. use_isotropic と inner Monte Carlo が存在しない。
4. 通常出力は CL, CL(avg), TR の 3 estimator に整理されている。
5. EI/CX は共通アルゴリズムとして実装されている。
6. EL は CL actual, CL(avg), TR で別々に評価される。
7. Sn/Sp/We/Wi が CSV 出力される。
8. step/block 統計から標準偏差と 95% 信頼区間が出る。
9. 粒子数、運動量、エネルギーの保存則確認 CSV が出る。
10. 断面積 CDF ファイルを差し替えて比較できる。
11. make OPENMP=1 で粒子ループを OpenMP 並列化できる。
12. make FAST=1 で bounds check なしの高速ビルドができる。
```
