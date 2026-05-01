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
  - pure-elastic では elastic-total と spin-exchange の `p(theta,E)` を runtime energy へ DCS-first で補間し、差を取ってから積分量と CDF を作る
- `Krstic/build_krstic_total_elastic_cdf.py`
  - `total elastic` 用の正規 builder
  - `sigma_t` と `sigma_mt`、I-kernel 系は Krstic integral fit を優先する
  - `scattering_angle` は DCS-first で作る: `p(theta,E)=2*pi*sin(theta)*d sigma/d Omega` を runtime energy へ補間し、それを積分して CDF / inverse-CDF にする
- `calc_I_kernel.py`
  - 比較用の補助スクリプト
  - Krstic pure DCS-first から `sigma_mt(E_cm)` を作り、I-kernel 系だけを再計算する
  - `cross_section`, `scattering_angle`, `angle_min` は入力 CDF から保持し、`reaction_rate`, `I_1_0`, `I_1_1_up`, `I_1_2_up2`, `sigv_max` だけを更新する
- `plot_cdf.py`
  - 既定では `Krstic/krstic_dd_total_elastic_integral_priority.cdf` を図示する
- `plot_cdf_test.py`
  - 既定ではルートの `dd_00_elastic_pure_el_angle_fixed.cdf` を図示する
- `debug_krstic_scattering_angle.py`
  - Bachmann baseline と Krstic total / pure の `scattering_angle` を同じ runtime energy 軸で比較する debug スクリプト
  - 代表 energy の線図、angle table から戻した transport ratio、per-energy diff サマリを出力する
- `check_energy_exchange_bias.py`
  - `dd_00_elastic.cdf` 互換 CDF の `reaction_rate`, `I_1_0`, `I_1_1*up`, `I_1_2*up^2` を読み、3d3v / 元コードの lookup 式で得られる弾性衝突のエネルギー交換率を中性粒子 Maxwell 分布で積分する検証スクリプト
  - 同温度・ゼロドリフトで `<sigma v DeltaE_n>` が 0 に近いかを見て、lookup table / 補間 / CDF 生成のバイアスを確認する
  - 符号は、正なら中性粒子側のエネルギー増加を表す。3d3v の plasma-side source `Wi` では反対符号になる
- `plot_energy_exchange_bias.py`
  - `check_energy_exchange_bias.py` の評価式を使い、`T_i = T_n` の温度掃引でイオン側エネルギー寄与 `-<sigma v DeltaE_n>` を図示する
  - 既定では 0.2--100 eV、温度クランプなしで、Krstic_total / Krstic_pure_el の2本を比較し、`energy_exchange_bias.txt` があれば実シミュレーションの `Wi` も重ねる
  - 図、理論曲線 CSV、シミュレーション換算 CSV を `figure/` に出力する

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
  - `sigma_t` と `sigma_mt`、I-kernel 系は Krstic integral fit を使う
  - `scattering_angle` は elastic-total DCS-first で再構成する
- `Krstic/krstic_dd_pure_dcs_compat.cdf`
  - `pure elastic` の正規データ
  - elastic-total / spin-exchange DCS を runtime energy へ DCS-first で補間し、`g_pure = g_el_total - g_se` から `cross_section`, `scattering_angle`, `reaction_rate`, `I_1_x`, `sigv_max` を再構成した full rebuild
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

## エネルギー交換率の検証

弾性衝突の track-length / lookup 評価では、実行時に散乱角を直接サンプリングせず、CDF 内の `I_1_x` テーブルを参照して運動量・エネルギー交換率を計算する。`check_energy_exchange_bias.py` はこの lookup 式を再現し、さらに中性粒子速度を Maxwell 分布で積分して、温度ペアごとの期待値

```text
<sigma v DeltaE_n>  [eV cm^3/s]
```

を評価する。たとえば 2.0 eV 同士のように中性粒子温度とイオン温度が等しく、ドリフトがない条件では、物理的には期待値が 0 に近いはずなので、残差は table 生成、補間、lookup 式のバイアス確認に使える。

```bash
# 検証用の既定: 温度クランプなし
python3 check_energy_exchange_bias.py --pairs 2,2 \
  --cdf Krstic/krstic_dd_total_elastic_integral_priority.cdf

# 元コードの TLMT_EL=0.9 eV/amu を再現する場合
python3 check_energy_exchange_bias.py --temperature-clamp --tlmt-el 0.9 \
  --pairs 1,1 2,2 --cdf dd_00_elastic.cdf
```

出力 CSV では、主結果の `lookup_gain_ev_cm3_s` に加えて、常に `lookup_gain_no_clamp_ev_cm3_s` と `lookup_gain_clamped_reference_ev_cm3_s` も出す。これにより、検証用の no-clamp 評価と、元コード再現用の clamped 評価を同じファイル内で比較できる。

図として確認する場合は、同じ期待値をイオン側の符号に反転してプロットする。

```bash
python3 plot_energy_exchange_bias.py
```

以前の4本比較を再生成したい場合は `--preset four` を使う。実シミュレーション点を重ねない純粋な lookup 曲線だけの比較にするなら、`--no-simulation-results` も指定する。

```bash
python3 plot_energy_exchange_bias.py --preset four --no-simulation-results \
  --output figure/energy_exchange_ion_equal_temperature_four_cdfs.png \
  --csv figure/energy_exchange_ion_equal_temperature_four_cdfs.csv
```

既定で参照する CDF は以下の2つである。Bachmann と Bachmann-Janev は既定の図には載せない。

- Krstic_total: `Krstic/krstic_dd_total_elastic_integral_priority.cdf`
- Krstic_pure_el: `Krstic/krstic_dd_pure_dcs_compat.cdf`

`--preset four` で参照する CDF は以下の4つである。

- Bachmann: `Bachmann/bachmann_dd_from_split_tables_compat.cdf`
- Bachmann-Janev: `dd_00_elastic_pure_el_angle_fixed.cdf`
- Krstic_total: `Krstic/krstic_dd_total_elastic_integral_priority.cdf`
- Krstic_pure_el: `Krstic/krstic_dd_pure_dcs_compat.cdf`

`energy_exchange_bias.txt` がある場合は、`Wi` と `stddev` を実シミュレーション結果として読む。`Wi` は `W/m^3` として扱い、既定では `n_i = 1.0e21 m^-3`, `n_n = 5.0e19 m^-3` を使って `eV cm^3/s` に換算して、点とエラーバーで重ねる。

既定の出力は以下である。

- `figure/energy_exchange_ion_equal_temperature.png`
- `figure/energy_exchange_ion_equal_temperature.csv`
- `figure/energy_exchange_simulation_converted.csv`

### `Etm` と `Tim`

CDF の `I_1_x(Etm, Tim)` は、テスト粒子である中性粒子の相対速度と、背景イオン Maxwell 分布の温度を引数にした事前積分テーブルである。元コードのコメントでは以下の定義になっている。

```text
Etm = 1/2 * mu * vt^2 = mu/mt * Et  [eV/amu]
Tim = mu/mi * Ti                    [eV/amu]
```

ここで `mu` は換算質量、`mt` は中性粒子質量、`mi` はイオン質量である。変数名 `tim` は時刻ではなく、この table temperature `Tim` を表す。D + D+ では `mu/mi = 1/2` なので、物理的なイオン温度 `Ti = 2.0 eV` は table 上では `tim = 1.0 eV/amu` になる。

元コードおよび 3d3v では、`I_1_0` と `I_1_1*up` の lookup だけに温度下限 `TLMT_EL = 0.9 eV/amu` をかけ、`I_1_2*up^2` と `reaction_rate` にはクランプ前の `tim` を使う。D + D+ の場合、これは物理温度では `Ti < 1.8 eV` のときだけ効く。したがって 2.0 eV 同士の検証では、このクランプは残差の原因にはならない。

## 注意点

- `calc_I_kernel.py` は full rebuild ではありません。比較用の I-kernel 差し替え器です。
- Krstic の DCS fit 入力は `E_cm` です。D + D では `E_cm = 0.5 * E_lab` を使います。
- `total elastic` の正式 CDF は、角度を DCS-first、積分量を integral-fit-first として分けます。DCS 由来の角度CDFから戻した transport ratio と integral fit の ratio が一致しない場合、その差は validation に残します。
- `pure elastic` の正式 CDF は、独立した pure の integral-fit が揃っていないため DCS-first full rebuild とします。elastic-total / spin-exchange の DCS を先に runtime energy へ補間し、差を取った pure kernel から積分量と角度 CDF を作ります。
- `g_pure = g_el_total - g_se` は、現行の手動入力 DCS 係数でも局所的に負になる点が残っています。これは raw OCR をそのまま使っているという意味ではなく、pure 再構成の整合性確認がまだ必要だという意味です。`build_krstic_angle_cdf.py` と `calc_I_kernel.py` は `negative-pure-policy` で扱いを切り替えます。
- より詳しい Krstic 側の変更履歴と生成物の説明は `Krstic/README.md` を参照してください。
