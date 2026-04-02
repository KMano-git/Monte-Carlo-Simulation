# 3D3V Lookup Check

`project/monte_simulations_ai/3d3v_lookup_check/` は、3D3V の非アナログ Monte Carlo 輸送に対して、
elastic source の score 方法を複数比較できるようにした作業ディレクトリです。

現在の主眼は、elastic について

- `A(EL)`
  - analog estimator
  - 実際に起きた EL 衝突の実現値 `-delta_E_el` を score
- `CL`
  - pure scatter-collision based collision estimator
  - アクセプトされた scatter collision ごとに、`CX/EL` の期待寄与を score
- `TR`
  - lookup-based track-length estimator
  - `I_1_0`, `I_1_1*up`, `I_1_2*up^2` から作る平均 rate を
    `pending_eff_time` に掛けて score

を mainline とし、elastic 比較用に

- `CLInner(EL)`
  - inner multi-sampling による elastic CL 比較値
- `TRInner(EL)`
  - inner multi-sampling による elastic TR 比較値
- `TRPretab(EL)`
  - offline Monte Carlo table の runtime lookup による elastic TR 比較値

を同時に出力します。

## Estimator の意味

### 1. `A(EL)`

elastic 実衝突時だけ score します。

```text
score += w * (-delta_E_el_actual)
```

### 2. `CL`

現在の輸送では、`CX` と `EL` の scatter collision が rejection を通過した瞬間が
collision score のタイミングです。

pure scatter-collision based CL として、各アクセプト scatter collision ごとに

```text
score_cl_cx += w * (R_cx / R_s) * s_cx
score_cl_el += w * (R_el / R_s) * s_el_bar
score_cl_ei += w * (R_a  / R_s) * s_a
```

を足します。ここで

- `R_s = R_cx + R_el`
- `s_cx`
  - その衝突点の sampled ion velocity に対する CX 寄与
- `s_el_bar`
  - elastic の 1 collision あたり期待寄与

です。

つまり、`CL` は actual channel を見て `CX` か `EL` のどちらかだけを足すのではなく、
アクセプトされた scatter collision ごとに両チャネルの期待寄与を score します。

### 3. `CLInner(EL)`

`CL` の elastic 項の比較用です。

- elastic 1 collision あたり期待寄与を inner multi-sampling で評価
- そのあと pure CL に合わせて `R_el / R_s` を掛ける

ので、`CL(EL)` と同じ scatter-collision based の定義で比較できます。

### 4. `TR`

lookup-based TR の mainline です。

elastic 項は `dd_00_elastic.cdf` の

- `reaction_rate`
- `I_1_0`
- `I_1_1*up`
- `I_1_2*up^2`

から unit-time expected source を作り、

```text
score_tr_el += s_c_el_lookup * pending_eff_time
```

を足します。

### 5. `TRInner(EL)`

elastic 項を inner multi-sampling で

```text
s_c_el = average_m [ n_i * sigma_el(E_rel_m) * v_rel_m * (-delta_E_el_m) ]
```

として評価し、`pending_eff_time` に掛ける比較用 TR です。

### 6. `TRPretab(EL)`

offline Monte Carlo で作った軽量 table を runtime lookup して、

```text
score_tr_pretab_el += s_c_el_table * pending_eff_time
```

を足す比較用 TR です。

## 実装の対応

### `code/scoring.f90`

主要 routine は次です。

- `score_analog_elastic(...)`
  - `A(EL)` を score
- `score_collision_estimator(...)`
  - pure scatter-collision based `CL`
- `score_inner_multi_collision_el(...)`
  - `CLInner(EL)`
- `score_track_length_estimator(...)`
  - lookup-based `TR`
- `score_inner_multi_track_length(...)`
  - `TRInner(EL)`
- `score_pretabulated_track_length(...)`
  - `TRPretab(EL)`

elastic の lookup kernel は

- `compute_el_lookup_coords(...)`
- `compute_el_lookup_rate_and_score(...)`

に分離しています。

### `code/main.f90`

実衝突時には

- `CX`
  - `CL`
  - `CLInner(EL)`
- `EL`
  - `CL`
  - `A(EL)`
  - `CLInner(EL)`

を score します。

flight flush 時には

- `TR`
- `TRInner(EL)`
- `TRPretab(EL)`

を score します。

### `code/io.f90`

`ntscrg.csv` とコンソール summary に、以下の比較結果を出します。

- `A_Q_el[W/m3]`
- `CL_Q_cx[W/m3]`
- `CL_Q_el[W/m3]`
- `CL_Q_ei[W/m3]`
- `CL_Q_total[W/m3]`
- `CLInner_Q_el[W/m3]`
- `TR_Q_cx[W/m3]`
- `TR_Q_el[W/m3]`
- `TR_Q_ei[W/m3]`
- `TR_Q_total[W/m3]`
- `TRInner_Q_el[W/m3]`
- `TRPretab_Q_el[W/m3]`

## 追加・変更した主なファイル

- `code/scoring.f90`
  - analog / pure CL / lookup TR と比較用 elastic scorer
- `code/main.f90`
  - scorer 呼び出し配線
- `code/io.f90`
  - CSV / summary 出力整理
- `code/data_types.f90`
  - score 項目の再編
- `code/tl_el_table_reader.f90`
  - pretabulated TR 用 table 読み書き
- `code/build_tl_el_table.f90`
  - offline table builder
- `Makefile`
  - main executable と table builder を build
- `run/input.nml`
  - runtime / builder 設定

## 設定

`run/input.nml` の `&simulation` では、少なくとも次を使います。

- `tl_el_inner_samples`
  - `CLInner(EL)` / `TRInner(EL)` のサンプル数
- `tl_el_table_file`
  - `TRPretab(EL)` 用 table ファイル名

table builder 側は `&tl_el_table_builder` を使います。

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

生成物:

- `run/monte_carlo_3d3v_natl`
- `run/build_tl_el_table`

### 2. pretab table 生成

```bash
cd /home/mano/project/monte_simulations_ai/3d3v_lookup_check/run
./build_tl_el_table
```

デフォルト出力:

- `tl_el_table.dat`
- `tl_el_table_stats.csv`

### 3. 本計算

```bash
cd /home/mano/project/monte_simulations_ai/3d3v_lookup_check/run
./monte_carlo_3d3v_natl
```

table が読めると起動時に

```text
TL EL pretab table loaded: tl_el_table.dat
```

と表示されます。

## 現時点の注意

- `CLInner(EL)`, `TRInner(EL)`, `TRPretab(EL)` の追加対象は elastic 項のみ
- `TR` mainline は `I_1_x` lookup ベース
- `TRPretab(EL)` は専用の軽量 table 形式で、`dd_00_elastic.cdf` そのものではない
- `A(EL)` は actual EL collision の analog score
- `CL` は pure scatter-collision based score
- `TR` は flight-based score

## このディレクトリで見たいこと

この構成にしてあることで、同じ条件下で

- analog EL の揺らぎがどの程度大きいか
- pure CL にするとどこまで分散が落ちるか
- lookup TR を mainline にしたときの振る舞い
- inner multi-sampling が lookup とどれくらい一致するか
- pretab TR が runtime lookup TR にどれくらい近いか

を切り分けやすくしています。
