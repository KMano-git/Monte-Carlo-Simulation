`el_data/` は、D + D+ の断面積テーブルまわりを扱う作業ディレクトリです。現在は「基準 CDF」「共通ユーティリティ」「Bachmann 系」「Krstic 系」「図出力」を分けて運用します。

## 置き場所の整理

- ルート直下
  - 基準 CDF と共通スクリプトを置く場所
  - `dd_00_elastic*.cdf` はここに残す
  - `figure/` はローカル図出力専用で、git 管理対象外
- `Bachmann/`
  - Bachmann 由来の抽出データと互換 CDF
- `Krstic/`
  - Krstic 由来の生成物と検証結果の正規置き場
  - Krstic 生成物はここを正本とし、`el_data/` 直下には重複配置しない
  - 係数の source of truth は `DpD_fit_memo_v2.md` の
    - Section 11: DCS 係数
    - Section 12: integral fit 係数
    - Section 2: integral reference table

## 主要ファイル

- `cdf_compat.py`
  - `dd_00_elastic.cdf` 互換の CDL テキストを読み書きする共通ユーティリティ
  - `compute_transport_tables(...)`
    - `scattering_angle` から `<1-cos(theta)> = R_theta` を作り、`sigma_tot * R_theta` を使う旧 surrogate 経路
  - `compute_transport_tables_from_sigma_momentum(...)`
    - `sigma_mt(E)` を直接与える新経路
- `build_elastic_cdf_compatible.py`
  - 分割した CSV/JSON 断面積データを `.cdf` 互換ファイルへ戻す汎用生成器
- `krstic_dcs.py`
  - Krstic の pure-elastic DCS evaluator と、`sigma_t / sigma_mt / sigma_vi / angle CDF` の再構成ユーティリティ
- `Krstic/generate_krstic_integral_data.py`
  - `DpD_fit_memo_v2.md` の Section 2 / Section 12 を読み、Krstic integral fit 系の JSON / CSV を再生成する
- `build_krstic_angle_cdf.py`
  - Krstic DCS から `cross_section`, `scattering_angle`, `reaction_rate`, `I_1_x`, `sigv_max` をまとめて再構成する full rebuild スクリプト
- `Krstic/build_krstic_total_elastic_cdf.py`
  - `total elastic` 用の正規 builder
  - `sigma_t` と `sigma_mt` は Krstic integral fit を優先し、`scattering_angle` は elastic-total DCS から再構成する
- `calc_I_kernel.py`
  - 比較用の補助スクリプト
  - Krstic DCS から `sigma_mt(E_cm)` を作り、I-kernel 系だけを再計算する
  - `cross_section`, `scattering_angle`, `angle_min` は入力 CDF から保持し、`reaction_rate`, `I_1_0`, `I_1_1_up`, `I_1_2_up2`, `sigv_max` だけを更新する
- `plot_cdf.py`
  - 既定では `Krstic/krstic_dd_pure_dcs_compat.cdf` を図示する
- `plot_cdf_test.py`
  - 既定ではルートの `dd_00_elastic_pure_el_angle_fixed.cdf` を図示する
- `debug_krstic_scattering_angle.py`
  - Bachmann baseline と Krstic total / pure の `scattering_angle` を同じ runtime energy 軸で比較する debug スクリプト
  - 代表 energy の線図、angle table から戻した transport ratio、per-energy diff サマリを出力する

## 基準 CDF の扱い

- `dd_00_elastic.cdf`
  - 現在の基準となる `dd_00_elastic.cdf` 互換データ
- `dd_00_elastic_pure_el.cdf`
- `dd_00_elastic_pure_el_angle.cdf`
- `dd_00_elastic_pure_el_angle_fixed.cdf`
- `dd_00_elastic_regen.cdf`
  - いずれもルート側の既存・派生 CDF 群
  - Krstic スクリプトの正規出力先ではない
  - とくに `dd_00_elastic_pure_el_angle_fixed.cdf` は非 Krstic の既存ファイルとして扱い、Krstic スクリプトで上書きしない

## Krstic 系ファイルの位置づけ

- `Krstic/krstic_dd_total_elastic_integral_priority.cdf`
  - `total elastic` の正規データ
  - `sigma_t` と `sigma_mt` は Krstic integral fit を使い、`scattering_angle` は elastic-total DCS から再構成する
- `Krstic/krstic_dd_pure_dcs_compat.cdf`
  - `pure elastic` の正規データ
  - pure-elastic DCS から `cross_section`, `scattering_angle`, `reaction_rate`, `I_1_x`, `sigv_max` を再構成した full rebuild
- `Krstic/krstic_dd_dcs_coeffs.json`
  - Section 11 の手動入力から作る共通 DCS 係数 JSON
  - total / pure の両 workflow で共有する

## いまの使い分け

- `total elastic` の正式データを使いたいとき
  - `Krstic/krstic_dd_total_elastic_integral_priority.cdf`
- `pure elastic` の正式データを使いたいとき
  - `Krstic/krstic_dd_pure_dcs_compat.cdf`
- 既存 CDF を下敷きにして I-kernel 差し替えの影響だけを見たいとき
  - `calc_I_kernel.py`
  - これは比較用であり、正式データの生成 workflow ではない

## 再生成

```bash
python3 Krstic/generate_krstic_integral_data.py
python3 Krstic/build_krstic_total_elastic_cdf.py
python3 build_krstic_angle_cdf.py
python3 debug_krstic_scattering_angle.py
python3 plot_cdf.py
python3 plot_cdf_test.py
```

比較用の hybrid 出力が必要な場合だけ、追加で以下を実行する:

```bash
python3 calc_I_kernel.py
```

## 注意点

- `calc_I_kernel.py` は full rebuild ではありません。比較用の I-kernel 差し替え器です。
- Krstic の DCS fit 入力は `E_cm` です。D + D では `E_cm = 0.5 * E_lab` を使います。
- `g_pure = g_el_total - g_se` は、現行の手動入力 DCS 係数でも局所的に負になる点が残っています。これは raw OCR をそのまま使っているという意味ではなく、pure 再構成の整合性確認がまだ必要だという意味です。`build_krstic_angle_cdf.py` と `calc_I_kernel.py` は `negative-pure-policy` で扱いを切り替えます。
- より詳しい Krstic 側の変更履歴と生成物の説明は `Krstic/README.md` を参照してください。
