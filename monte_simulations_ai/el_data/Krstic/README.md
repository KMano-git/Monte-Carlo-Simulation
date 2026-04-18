Krstic D+ + D working files

このディレクトリは、Krstic 由来の生成物と検証結果の正規置き場です。`el_data/` 直下に一時的に置かれていた Krstic 生成物は、現在はここに集約します。

## 正規 workflow

1. `generate_krstic_integral_data.py`
   - `../DpD_fit_memo_v2.md` の Section 2 / Section 12 を読み、
     integral fit / reference table を JSON / CSV に整形する
2. `build_krstic_total_elastic_cdf.py`
   - `total elastic` の正規 workflow
   - `sigma_t` / `sigma_mt` は integral fit を優先
   - `scattering_angle` は elastic-total DCS から再構成
3. `../build_krstic_angle_cdf.py`
   - `pure elastic` の正規 workflow
   - `sigma_t`, `sigma_mt`, `scattering_angle`, `reaction_rate`, `I_1_x`, `sigv_max` を pure-elastic DCS から再構成する
4. `../calc_I_kernel.py`
   - 比較用の補助 workflow
   - 入力 CDF の `cross_section`, `scattering_angle`, `angle_min` を保持したまま、I-kernel 系だけを Krstic `sigma_mt` で差し替える
   - 正式データの生成には使わない

## 主要出力

- `krstic_dd_total_elastic_integral_priority.cdf`
  - `total elastic` の正規データ
  - `sigma_t` / `sigma_mt` は integral fit 優先
  - `scattering_angle` は elastic-total DCS
- `krstic_dd_pure_dcs_compat.cdf`
  - `pure elastic` の正規データ
  - pure-elastic DCS から `cross_section`, `scattering_angle`, `reaction_rate`, `I_1_x`, `sigv_max` を再構成した full rebuild
- `krstic_dd_dcs_coeffs.json`
  - Section 11 から作る共通 DCS 係数 JSON
  - total / pure の両 workflow で共有する

## 付随ファイル

- `krstic_dd_integral_coeffs.json`
- `krstic_dd_integral_reference_points.csv`
- `krstic_dd_integral_cdf_grid_101.csv`
- `krstic_dd_integral_angle_grid_51.csv`
- `krstic_total_elastic_validation.csv`
- `krstic_total_elastic_validation.json`
- `krstic_total_elastic_scattering_angle_compat.csv`
- `krstic_total_elastic_transport.csv`
- `krstic_dd_dcs_coeffs.json`
- `krstic_pure_dcs_validation.csv`
- `krstic_pure_dcs_validation.json`
- `krstic_pure_scattering_angle_compat.csv`
- `krstic_pure_transport_from_dcs.csv`

## スクリプトの役割

- `generate_krstic_integral_data.py`
  - Section 2 の reference table と Section 12 の manual fitting parameters を整理
- `build_krstic_total_elastic_cdf.py`
  - total-elastic 正式版を再生成
- `../krstic_dcs.py`
  - Krstic DCS evaluator
- `../build_krstic_angle_cdf.py`
  - pure-elastic DCS から full rebuild を作る
- `../calc_I_kernel.py`
  - pure-elastic DCS から `sigma_mt` を作り、I-kernel 系だけを差し替える比較用スクリプト

## 再生成

```bash
python3 generate_krstic_integral_data.py
python3 build_krstic_total_elastic_cdf.py
python3 ../build_krstic_angle_cdf.py
```

比較用の hybrid 出力を確認したい場合だけ、追加で以下を実行する:

```bash
python3 ../calc_I_kernel.py
```

## 現在の注意点

- `krstic_dd_total_elastic_integral_priority.cdf` と `krstic_dd_pure_dcs_compat.cdf` が、現在の 2 本の正式データです。
- `../calc_I_kernel.py` は比較用であり、正式な pure-elastic データ生成 workflow ではありません。
- integral fit の source of truth は `../DpD_fit_memo_v2.md` の Section 12、DCS の source of truth は Section 11 です。
- pure-elastic DCS は `g_pure(theta, E) = g_el_total(theta, E) - g_se(theta, E)` で作っています。
- 現行の手動入力 DCS 係数でも、一部エネルギーで `g_pure` が局所的に負になります。これは raw OCR をそのまま使っているという意味ではなく、pure 再構成の整合性確認がまだ必要だという意味です。`negative-pure-policy=warn-clip` は実装・比較用には便利ですが、最終 production 用データとして扱う前に見直しが必要です。
