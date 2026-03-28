# VSCode の ChatGPT / Codex に渡す最小設計書

## 目的

この作業の主目的は、既存の `scoring.f90` に対して、elastic 項の TL 評価方法として次の 2 方式を実装・比較できるようにすることである。

1. **inner multi-sampling TL**

   * TL の内部で 1 回だけ仮想衝突を起こす代わりに、複数回サンプリングして平均寄与を評価する方式

2. **pre-tabulated TL**

   * 事前に Monte Carlo サンプリングで平均寄与を計算し、その結果を table に保存して、実行時は lookup する方式

ここでは大規模な全面改修は行わず、**現行コードを温存しつつ、2 方式を追加して比較可能にする**ことを目標とする。

---

## 現状認識

`scoring.f90` の現行実装では、

* `score_collision_estimator` は実衝突時に energy transfer を加算する
* `score_track_length_estimator` は `pending_eff_time` に rate を掛けて加算する
* 特に EL 項では、内部で 1 回の仮想衝突をサンプルして `delta_E_el_dummy` を作る
* `score_table_lookup_track_length` は EL 項を `I_1_0`, `I_1_1`, `I_1_2` table から直接計算する

という構造である。つまり今の naive TL は、TL の外枠の中で EL 寄与を 1 サンプルの仮想 collision で見積もっている。 fileciteturn8file0

---

## 今回やること

今回やるのは次の 3 点に限定する。

### A. 現行 naive TL をベースに、inner multi-sampling 版を追加する

現行の elastic 寄与評価

* 1 回だけイオン速度サンプリング
* 1 回だけ散乱角サンプリング
* 1 回だけ `delta_E_el_dummy`

を、`N_inner` 回繰り返して平均化する版を追加する。

目標式:

```text
s_c_el = average over m=1..N_inner of [ n_i * sigma_el(E_rel_m) * v_rel_m * (-delta_E_el_m) ]
score_tl_el += s_c_el * pending_eff_time
```

これは、今の naive TL を最小限の変更で「複数回仮想衝突平均版」にしたものである。

---

### B. 事前サンプリング table を作る

elastic 項について、実行前に offline Monte Carlo で平均寄与を計算し、table に保存する。

最初の段階では、複雑にしすぎず、lookup 変数は現行 lookup にできるだけ合わせる。

第一候補:

* `etm`
* `tim`

必要なら後で

* `ut`
* `T_i`
* `E_kin`

を加える。

各 grid 点で多数回サンプリングして、たとえば

```text
mean_el_rate(etm, tim) = < n_i * sigma_el * v_rel * (-delta_E_el) >
```

を保存する。

---

### C. 実行時に 3 方式を比較できるようにする

比較対象は最小限、次の 3 つでよい。

1. 現行 TL naive
2. inner multi-sampling TL
3. pre-tabulated TL

余裕があれば基準として CL も残す。

---

# VSCode の ChatGPT / Codex に渡す最初の指示

```text
Fortran コードの最小差分改修を手伝ってください。
目的は scoring.f90 の elastic TL 項について、

1. 現行 naive TL
2. inner multi-sampling TL
3. pre-tabulated TL

を比較できるようにすることです。

重要:
- 大規模設計変更はしない
- 既存コードは壊さない
- 現行 naive 実装は残す
- 新方式は追加実装にする
- まずは elastic 項だけを対象にする
- 変更は小さく、各ステップで diff を出す

最初にやること:
1. scoring.f90 を読んで、elastic TL 項だけのデータフローを整理
2. delta_E_el_dummy が作られる処理を関数化できるか確認
3. inner multi-sampling 版を最小変更で入れる案を出す
4. offline table 作成ツールの最小設計を出す
```

---

# 実装ステップ

## Step 1: 現行 EL naive の関数化

まず、現行 `score_track_length_estimator` の EL 項を切り出す。

目標:

```fortran
real(dp) function tl_el_sample_once(...)
```

この関数は 1 回分の

* ion velocity sampling
* scattering angle sampling
* `delta_E_el_dummy`
* `n_i * sigma_el * v_rel * (-delta_E_el_dummy)`

を返す。

この段階では**数値挙動を変えない**。

---

## Step 2: inner multi-sampling 版の追加

次に、上の `tl_el_sample_once(...)` を `N_inner` 回呼び出して平均する関数を追加する。

目標:

```fortran
real(dp) function tl_el_sample_avg(..., n_inner)
```

実装イメージ:

```fortran
sum_rate = 0d0
Do m = 1, n_inner
   sum_rate = sum_rate + tl_el_sample_once(...)
End Do
tl_el_sample_avg = sum_rate / real(n_inner, dp)
```

そして新しい scoring subroutine では

```fortran
s_c_el = tl_el_sample_avg(..., n_inner)
score%tl_el = score%tl_el + s_c_el * pending_eff_time
```

とする。

---

## Step 3: offline table 作成ツールの追加

新規ファイルを 1 つ追加する。本プログラムの前に実行する別プログラムとして設計する。
統計精度を確保するために十分な回数の平均を取ること。

候補名:

* `build_tl_el_table.f90`

役割:

* `(etm, tim)` grid をループ
* 各点で `N_table_samples` 回サンプリング
* 平均 elastic rate を計算するのに使用
* 現行のCLによる結果を用いる
* dd_00_elastic.cdfのフォーマットに則って書き出す

## Step 4: table reader の追加

既存の `cdf_reader.f90` に無理に混ぜる必要はないが、基本的には同じ構造を用いたほうが楽だと推測される。

候補名:

* `tl_el_table_reader.f90`

機能:

* テーブルを書き出したファイルの読み込み
* nearest / bilinear interpolation
* `get_tl_el_rate_from_table(etm, tim)` を返す

---

## Step 5: pre-tabulated TL subroutine の追加

`score_table_lookup_track_length` をいじりすぎず、別の subroutine を追加する。

候補名:

* `score_pretabulated_track_length`

elastic 項だけ table を使い、EI/CX は現行のままでもよい。

実装イメージ:

```fortran
s_c_el = get_tl_el_rate_from_table(etm, tim)
score%tl_pretab_el = score%tl_pretab_el + s_c_el * pending_eff_time
```

---

# まず比較するもの

最初の比較はこれだけでよい。

### Case 1: 現行 naive

* 1 回仮想衝突

### Case 2: inner multi-sampling

* 10回、50回、100回など

### Case 3: pre-tabulated

* offline で十分大きいサンプル数で平均化した table を利用

見る量:

* EL score の平均
* EL score の分散
* 実行時間

# これだけは確認すること

1. `p%weight` が TL 側でどこで掛かっているか
2. `pending_eff_time` がどこで作られ、どの単位か
3. inner multi-sampling で RNG の使い方をどうするか
4. table の独立変数を `(etm, tim)` で十分とみなせるか

---

# 今回の方針の一言要約

今回の実装は、

* **今の naive TL を少しだけ拡張して複数回平均版を作る**
* **同じ平均を offline 化して table にする**
* **両者を比較する**

という最小構成で進める。

大規模な estimator 再設計ではなく、**EL 寄与評価の方法だけを比較対象として差し替える**のが今回の主眼である。
