Krstic D+ + D working files

このディレクトリは、Krstic 由来の生成物と検証結果の正規置き場です。`el_data/` 直下に一時的に置かれていた Krstic 生成物は、現在はここに集約します。

## 正規 workflow

1. `generate_krstic_integral_data.py`
   - `../DpD_fit_memo_v2.md` の Section 2 / Section 12 を読み、
     integral fit / reference table を JSON / CSV に整形する
2. `build_krstic_total_elastic_cdf.py`
   - `total elastic` の正規 workflow
   - `sigma_t` / `sigma_mt` と I-kernel 系は integral fit を優先
   - `scattering_angle` は elastic-total DCS を runtime energy へ DCS-first で補間して再構成
   - DCS-first は `p(theta,E)=2*pi*sin(theta)*d sigma/d Omega` を log(E)-log(p) 補間し、その後で CDF / inverse-CDF を作る方針
3. `../build_krstic_angle_cdf.py`
   - `pure elastic` の正規 workflow
   - `sigma_t`, `sigma_mt`, `scattering_angle`, `reaction_rate`, `I_1_x`, `sigv_max` を pure-elastic DCS から再構成する
   - elastic-total / spin-exchange の `p(theta,E)` を runtime energy へ DCS-first で補間し、差を取った pure kernel から積分量と CDF を作る
4. `../calc_I_kernel.py`
   - 比較用の補助 workflow
   - 入力 CDF の `cross_section`, `scattering_angle`, `angle_min` を保持したまま、I-kernel 系だけを Krstic pure DCS-first `sigma_mt` で差し替える
   - 正式データの生成には使わない
5. `../debug_krstic_scattering_angle.py`
   - Bachmann baseline と Krstic total / pure の angle table を比較する debug workflow
   - representative line plots と per-energy diff summary を出力する

## 主要出力

- `krstic_dd_total_elastic_integral_priority.cdf`
  - `total elastic` の正規データ
  - `sigma_t` / `sigma_mt` と I-kernel 系は integral fit 優先
  - `scattering_angle` は elastic-total DCS-first
- `krstic_dd_pure_dcs_compat.cdf`
  - `pure elastic` の正規データ
  - elastic-total / spin-exchange DCS を runtime energy へ DCS-first で補間し、`g_pure = g_el_total - g_se` から `cross_section`, `scattering_angle`, `reaction_rate`, `I_1_x`, `sigv_max` を再構成した full rebuild
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
  - 積分量は integral-fit-first、角度分布は DCS-first として扱う
- `../krstic_dcs.py`
  - Krstic DCS evaluator
- `../build_krstic_angle_cdf.py`
  - pure-elastic DCS から full rebuild を作る
- `../calc_I_kernel.py`
  - pure-elastic DCS-first から `sigma_mt` を作り、I-kernel 系だけを差し替える比較用スクリプト
- `../debug_krstic_scattering_angle.py`
  - Bachmann baseline と Krstic angle table の差を数値化する

## 再生成

```bash
python3 generate_krstic_integral_data.py
python3 build_krstic_total_elastic_cdf.py
python3 ../build_krstic_angle_cdf.py
python3 ../debug_krstic_scattering_angle.py
```

比較用の hybrid 出力を確認したい場合だけ、追加で以下を実行する:

```bash
python3 ../calc_I_kernel.py
```

## 現在の注意点

- `krstic_dd_total_elastic_integral_priority.cdf` と `krstic_dd_pure_dcs_compat.cdf` が、現在の 2 本の正式データです。
- `../calc_I_kernel.py` は比較用であり、正式な pure-elastic データ生成 workflow ではありません。
- `../debug_krstic_scattering_angle.py` は、Bachmann との差を angle table レベルで追うための debug 補助です。
- integral fit の source of truth は `../DpD_fit_memo_v2.md` の Section 12、DCS の source of truth は Section 11 です。
- total-elastic の正式データでは、角度CDFを integral fit に無理に合わせません。角度は DCS-first、`sigma_t` / `sigma_mt` / I-kernel 系は integral-fit-first として validation に差を残します。
- pure-elastic の正式データは、独立した pure の integral-fit が揃っていないため DCS-first full rebuild とします。elastic-total / spin-exchange の DCS を先に runtime energy へ補間し、`g_pure(theta, E) = g_el_total(theta, E) - g_se(theta, E)` から積分量と角度 CDF を作ります。
- 現行の手動入力 DCS 係数でも、一部エネルギーで `g_pure` が局所的に負になります。これは raw OCR をそのまま使っているという意味ではなく、pure 再構成の整合性確認がまだ必要だという意味です。`negative-pure-policy=warn-clip` は実装・比較用には便利ですが、最終 production 用データとして扱う前に見直しが必要です。
