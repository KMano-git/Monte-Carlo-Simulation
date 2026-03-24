Event-driven Monte Carlo Simulation

愚直な track-length estimator とテーブルルックアップの track-length を比較したい

## 判明した問題点メモ

### 1. 修正済みの実装上の問題

- Track-Length scorer が粒子本体の RNG を直接進めていた。
  これにより、diagnostics を有効化すると実際の衝突列まで変わっていた。
  現在は scorer 内で local copy の RNG を使うよう修正済み。

- Track-Length estimator の呼び出しタイミングが「区間ごと」になっていた。
  現在は「実衝突確定時に pending 区間を flush、最後だけ終了時 flush」に変更済み。

- `t_col == time_remaining` のとき、ステップ末尾の衝突評価が落ちる問題があった。
  `time_remaining > 0` 依存を外して修正済み。

- `enable_cx` または `enable_el` の片側だけを無効化した場合、
  受理済み衝突が `COLL_NONE` になる経路があった。
  現在は enable されている断面積だけで acceptance/type selection を行うよう修正済み。

- EL の CDF 入力軸は `specific_energy [eV/amu]` だが、
  実装側では `E_rel [eV]` をそのまま入れていた。
  現在は断面積・散乱角 lookup ともに `0.5 * E_rel` を使うよう修正済み。

- `cdf_reader.f90` で `xs_mult` をハードコードしていた。
  現在は CDF から `xs_mult` を読み出すよう修正済み。

### 2. 単位変換そのものは主因ではなかった

- `1e-6 / 1e-8 / 1e-10` の CGS -> MKS 変換は、
  CDF の `xs_mult` および DEGAS 側の `xwmlt` と整合している。

- `calc_I_kernel.py` 内の `c0 = 1.3891494e6 cm/s` についても、
  `dd_00_elastic_pure_el_angle_fixed.cdf` を再積分すると
  `reaction_rate`, `I_1_0`, `I_1_1*up`, `I_1_2*up^2` が機械精度レベルで一致した。

- したがって、単純な unit conversion の誤りが
  `TL Lookup EL` の大差の主因ではない。

### 3. `dd_00_elastic.cdf` の内部不整合

- `dd_00_elastic.cdf` について、
  保存済みの `cross_section` と `scattering_angle` から
  `calc_I_kernel.py` と同じ手順で `reaction_rate / I_1_x` を再積分すると、
  `reaction_rate` はほぼ再現できる一方、`I_1_x` は再現できない。

- 再積分との比較結果:
  - `reaction_rate`: mean relative error 約 `3.0e-3`
  - `I_1_0`: mean relative error 約 `4.65e-1`
  - `I_1_1*up`: mean relative error 約 `4.63e-1`
  - `I_1_2*up^2`: mean relative error 約 `4.51e-1`

- 代表点では、再積分値 / 保存値 がだいたい `0.42 - 0.48` 程度であり、
  `dd_00_elastic.cdf` の保存済み `I_1_x` は、
  同じ CDF 内の `sigma_tot + theta(R,E)` からは出てこない。

### 4. EL の大差は 2 段構造

- `dd_00_elastic.cdf` をそのまま使った EL-only 短時間再現では、
  `TL_EL ~= 1.48e6 W/m3`, `CL_EL ~= 1.17e6 W/m3` に対して
  `TLLookup_EL ~= 3.61e7 W/m3` となった。

- `dd_00_elastic.cdf` の `sigma_tot + theta` から
  `I_1_x` だけ再生成した一時 CDF に差し替えると、
  `TLLookup_EL ~= 1.25e7 W/m3` まで低下した。

- つまり、元の `dd_00_elastic.cdf` に含まれる `I_1_x` の不整合だけで
  約 `2.9x` の過大化が入っている。

- ただし、それでも `TL_EL` / `CL_EL` よりなお大きく、
  残差はさらに別要因を含む。

### 5. 論文側の注意点と整合する残差

- Bachmann & Reiter (1995) では、
  cumulative-angle 的な近似は低エネルギー側で
  momentum / energy transfer cross section に
  factor `4 - 5` 級のずれを生じうると記述されている。

- 本コードでも、`fixed` 側のように
  `sigma_tot + theta -> I_1_x` が内部整合している CDF を用いても、
  EL については
  `TL/CL` と `TLLookup` がまだ大きくは一致しない。

- よって現在の理解では、
  EL の大差は次の 2 つが重なっている:
  1. `dd_00_elastic.cdf` 内の保存済み `I_1_x` が `sigma_tot + theta` と整合していない
  2. angle-sampling に基づく愚直 TL/CL と、`I_1_x` に基づく lookup EL は、
     低エネルギー域ではそもそも数値的に大きくずれうる

## 現時点の結論

- 実装上の obvious bug は一通り修正済み。
- 単位変換の誤りは主因ではない。
- `dd_00_elastic.cdf` を基準にする場合、
  まず「保存済み `I_1_x` を信用してよいか」を分けて考える必要がある。
- 今後の検証では、
  `saved I_1_x` と `sigma_tot + theta から再生成した I_1_x` を
  切り替えて比較できるようにしておくのが有効。
