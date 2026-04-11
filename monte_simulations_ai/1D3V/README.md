# 1D3V Monte Carlo Simulation

`1D3V` は、空間の扱いを `x` のみに限定しつつ、粒子の位置と速度は内部的に `x,y,z,vx,vy,vz` を追跡する中性粒子輸送コードです。実装は `3d3v_lookup` 系の event-driven / null-collision / non-analog 重み更新を基準にそろえています。

## 特徴

- `x` 方向だけを空間ビニングし、`y,z` は自由飛行の内部座標として保持
- CX / EL / EI を考慮
- 粒子固有 RNG による再現性確保
- EL CDF 軸を DEGAS と整合する `E_rel` に修正し、`xs_mult` 読み出しを反映
- 標準では `CL` と `TL` を出力し、`enable_tl_lookup=.true.` のときだけ `TLLu` を追加

## ディレクトリ構成

```text
1D3V/
├── code/
├── module/
├── run/
├── Makefile
├── README.md
└── plot_results.py
```

## ビルドと実行

```bash
make
cd run
./monte_carlo_1d3v
cd ..
python3 plot_results.py
```

OpenMP フラグ付きビルド:

```bash
make OPENMP=1
```

## 入力ファイル

`run/input.nml` は次の 4 セクションを持ちます。

- `simulation`
  - `n_particles, n_steps, n_x_bins, dt`
  - `x_min, x_max, gamma_in`
  - `enable_cx, enable_el, enable_ei, use_isotropic, enable_tl_lookup`
  - `weight_min`
  - `cdf_file, output_ntscrg, output_profile, output_hist, output_deltaE_hist`
- `plasma_nml`
  - `n_i, T_i, n_e, T_e, u_x, u_y, u_z`
- `particle_init`
  - `x_init, y_init, z_init, E_init, T_init, init_mode`
- `diagnostics`
  - `output_interval`
  - `n_hist_bins, E_hist_min, E_hist_max, hist_timing`
  - `n_dE_bins, dE_hist_min, dE_hist_max, dE_collect_steps`

## 出力ファイル

- `run/ntscrg.csv`
  - 1 ステップ 1 行の全体時系列
  - `CL_Q_*[W/m3]`, `TL_Q_*[W/m3]`
  - `enable_tl_lookup=.true.` のときだけ `TLLu_Q_*[W/m3]`
- `run/profile_x.csv`
  - `x` ビンごとの最終空間分布
  - `neutral_density[m-3]`
  - `CL_Q_*[W/m3]`, `TL_Q_*[W/m3]`
  - `collision_count_*_per_history[-]`
  - `enable_tl_lookup=.true.` のときだけ `TLLu_Q_*[W/m3]`
- `run/energy_hist.csv`
  - 生存粒子エネルギー分布の履歴
- `run/deltaE_hist.csv`
  - `enable_tl_lookup=.true.` のときのみ出力

## 正規化

- 時系列の全体ソース:
  - `Q_global = score_step[J] * gamma_in / (N_particles * dt * Lx)`
- 空間プロファイルの局所パワー密度:
  - `Q_x = score_bin[J] * gamma_in / (N_particles * total_time * dx)`
- 中性粒子密度:
  - `n_x = residence_eff[s] * gamma_in / (N_particles * dx)`

ここで `residence_eff` は non-analog 重みを含んだ有効滞在時間です。

## Lookup モードについて

`dd_00_elastic.cdf` に含まれる `I_1_x` には内部不整合の可能性があるため、`TLLu` は標準では無効です。比較診断をしたいときだけ `enable_tl_lookup=.true.` を使ってください。

## 可視化

`plot_results.py` は以下を作成します。

- `run/figure/1d3v_summary.png`
- `run/figure/1d3v_profiles.png`

`ntscrg.csv` と `profile_x.csv` の列を自動判別し、lookup 列があれば同じ図に重ねて描画します。
