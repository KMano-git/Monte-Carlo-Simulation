# 3D3V Lookup Check

`scoring.f90` の elastic TL 項について、既存実装を壊さずに複数方式を比較できるようにした作業ディレクトリです。

現在このディレクトリでは、以下の EL 評価を同時に比較できます。

- `CL`
  - 実衝突時の collision estimator
- `TL`
  - 既存の naive track-length estimator
  - TL 内で 1 回だけ仮想 elastic collision をサンプル
- `TLInner(EL)`
  - inner multi-sampling 版
  - TL 内で複数回サンプリングして平均
- `TLPretab(EL)`
  - offline Monte Carlo で作成した table を runtime lookup
- `TLLookup`
  - 既存の `I_1_x` ベース lookup 実装

## この改修でやったこと

### 1. task.md の内容を精査して方針を補正

`task.md` の主旨はそのまま採用しつつ、実装上は次の 2 点を補正しています。

- `pre-tabulated TL` は既存の `I_1_x` ベース `TLLookup` とは別物として追加
  - 既存 lookup はそのまま残し、比較対象として併存
- offline table は `CL` の平均ではなく、
  `naive TL` の EL 1 サンプル核
  `n_i * sigma_el * v_rel * (-delta_E_el_dummy)`
  の平均として実装
  - これにより `TL naive` と `TL inner` の延長として解釈しやすくしています

また、TL 側の重みの扱いは既存コードに合わせています。

- `p%weight` は scorer 内で直接掛けず、
  `pending_eff_time` にすでに織り込み済み
- `pending_eff_time` は `compute_effective_track_time()` で生成
  - EI loss がある場合は
    `weight / loss_rate * (1 - exp(-loss_rate * dt))`
  - loss が無視できる場合は `weight * dt`

## 実装内容

### EL naive TL の関数化

`code/scoring.f90` で、既存 naive TL の EL 部分を再利用できるように分離しました。

- `tl_el_rate_from_ion_sample(...)`
  - 背景イオン速度が与えられたときの EL 1 サンプル rate を返す
- `tl_el_sample_once(...)`
  - naive TL と同じ 1 回サンプル
- `tl_el_sample_avg(..., n_inner)`
  - inner multi-sampling の平均

既存の `score_track_length_estimator(...)` は、元の 1 サンプル挙動を保ったまま上記関数を使う形にしています。

### inner multi-sampling TL の追加

`code/scoring.f90` に以下を追加しました。

- `score_inner_multi_track_length(...)`

この scorer は elastic 項だけを評価し、

```text
s_c_el = average_m [ n_i * sigma_el(E_rel_m) * v_rel_m * (-delta_E_el_m) ]
score += s_c_el * pending_eff_time
```

に対応します。

### pre-tabulated TL の追加

専用の軽量 table 形式を追加しました。

- `code/tl_el_table_reader.f90`
  - table の read/write
  - default grid 生成
  - bilinear interpolation
- `code/build_tl_el_table.f90`
  - offline Monte Carlo で EL table を生成
- `code/scoring.f90`
  - `score_pretabulated_track_length(...)`

lookup 変数は最初の task 方針に合わせて、まずは `(etm, tim)` を使っています。

- `etm`
  - drift ベースの相対エネルギー座標
- `tim`
  - 背景イオン温度座標

### main への配線

`code/main.f90` では、pending TL flush のたびに以下を実行します。

1. `TL` naive
2. `TLInner(EL)`
3. `TLPretab(EL)`
4. 既存 `TLLookup`

これにより、同一の粒子履歴に対して各 estimator を同時比較できます。

### 出力の追加

`ntscrg.csv` に以下の列を追加しました。

- `TLInner_Q_el[W/m3]`
- `TLPretab_Q_el[W/m3]`

コンソール最終出力にも以下を追加しています。

- `TLInner(EL)`
- `TLPretab(EL)`

## 追加・変更した主なファイル

- `code/scoring.f90`
  - EL naive の関数化
  - inner multi-sampling scorer 追加
  - pre-tabulated scorer 追加
- `code/tl_el_table_reader.f90`
  - 新規
- `code/build_tl_el_table.f90`
  - 新規
- `code/main.f90`
  - table 読込
  - 新 scorer 呼び出し
- `code/io.f90`
  - namelist 拡張
  - CSV / 最終表示拡張
- `code/data_types.f90`
  - 新設定・新スコア項目追加
- `Makefile`
  - main executable に加えて table builder も build
- `run/input.nml`
  - 新パラメータ追加

## 新しく追加した設定

`run/input.nml` の `&simulation` に以下を追加しました。

- `tl_el_inner_samples`
  - inner multi-sampling のサンプル数
- `tl_el_table_file`
  - pre-tabulated TL 用 table ファイル名

また、table builder 用に新しい namelist を追加しています。

- `&tl_el_table_builder`
  - `n_table_samples`
  - `output_file`
  - `stats_output_file`
  - `rel_stderr_thresholds`

## ビルドと実行

### 1. ビルド

```bash
cd /home/mano/project/monte_simulations_ai/3d3v_lookup_check
make all
```

生成される実行ファイル:

- `run/monte_carlo_3d3v_natl`
- `run/build_tl_el_table`

### 2. pre-tabulated TL 用 table 生成

```bash
cd /home/mano/project/monte_simulations_ai/3d3v_lookup_check/run
./build_tl_el_table
```

デフォルトでは以下を出力します。

- `tl_el_table.dat`
  - runtime lookup 用の平均 table
- `tl_el_table_stats.csv`
  - 各グリッド点の統計情報
  - `mean_rate`, `stddev`, `variance`, `stderr`, `rel_stderr`

`rel_stderr = stderr / abs(mean)` です。
ただし、近平衡のように平均値そのものが 0 に近い点では `rel_stderr` が非常に大きくなりやすいので、
その場合は `stderr` の絶対値も併せて見るのが安全です。

builder の最後には、`input.nml` の `rel_stderr_thresholds` に基づいて、たとえば次のような集計も表示されます。

```text
rel stderr >  1.00%:  2601 / 2601
rel stderr >  5.00%:  1488 / 2601
rel stderr > 10.00%:   226 / 2601
max rel stderr: 4.8657E+01
max stderr:     7.2952E-11
```

### 3. 本計算の実行

```bash
cd /home/mano/project/monte_simulations_ai/3d3v_lookup_check/run
./monte_carlo_3d3v_natl
```

起動時に table が読めた場合は、

```text
TL EL pretab table loaded: tl_el_table.dat
```

と表示されます。

## 現時点の制約

- `inner multi-sampling` と `pre-tabulated` の追加対象は elastic 項のみ
- `TLPretab(EL)` は専用の軽量 table 形式を使用
  - 既存 `dd_00_elastic.cdf` の完全互換形式にはしていません
- `TLLookup` の既存 `I_1_x` ベース経路はそのまま残してあり、
  今回追加した `TLPretab(EL)` とは物理的な意味づけが異なります
- 分散や統計誤差の自動集計は未実装
  - simulation 本体側では未実装
  - ただし table builder 側では `tl_el_table_stats.csv` に
    `stderr` と `rel_stderr` を出力するようにしています

## 今回の実装の狙い

大規模な設計変更はせず、既存コードを温存したまま、

- 既存 naive TL
- inner multi-sampling TL
- pre-tabulated TL
- 既存 `I_1_x` lookup

を同じ条件で見比べられる状態にすることを優先しています。

そのため、今後の検証では以下を切り分けやすくなっています。

- naive TL の 1 サンプルノイズがどの程度支配的か
- inner averaging でどこまで安定化するか
- pre-tabulated にしたときに naive TL 平均へ近づくか
- 既存 `TLLookup` の差が estimator ノイズではなく kernel 側にあるか
