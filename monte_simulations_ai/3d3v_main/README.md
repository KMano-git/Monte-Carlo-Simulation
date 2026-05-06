# 3D3V Monte Carlo Simulation

`project/monte_simulations_ai/3d3v_main/`はこれまでの非アナログモンテカルロシミュレーションをまとめたコードである。
外部発表のため、Neut2Dのスラブコードとして完成させるためのディレクトリである。

元にしたコードは`project/monte_simulations_ai/3d3v_lookup_check/`である。

## 変更点

1. 旧Analog estimator $\rightarrow$ 新Collision estimator (名称の統一)
2. 旧Collision estimator $\rightarrow$ 新Collision estimator (averaging)
3. Track-length (pretable) はflag 変数 compare_estimator = .true.でのみ使用
4. 最終的には3つのCL,CL(avg),TRの3つのestimatorの出力を見れれば良い
5. Inner monte carlo は削除
6. 推定器の推定量を増やす：運動量移行率、粒子ソースなど
7. データの後処理を組み込んでみる
   1. 各推定結果の標準偏差を出したい
   2. 各保存則の収支
      1. 粒子数
      2. 運動量
      3. エネルギー
   3. 各推定器の信頼区間とやらもあったらいいかも
8. 等方散乱フラグは完全に使用しないので消す（その処理も同様）

これらのコードの当面の目的は、使用する断面積データによる違いを観察することである。

## 設計書

実装方針、推定器の定義、EI/CX/EL の source 計算式、出力仕様は [DESIGN.md](DESIGN.md) にまとめる。

## CX モデルの比較

`run/input.nml` の `cx_model` で CX の扱いを切り替えられる。

```fortran
cx_model = 0  ! 現行 3D3V sampled ion 版
cx_model = 1  ! src_read/monte 準拠版
```

`cx_model = 0` は、衝突判定時に背景イオン速度 `vi` を Maxwell 分布から sample し、

```text
v_rel = |v0 - vi|
E_rel = 0.25 * m_D * v_rel^2 / e
R_CX  = n_i * sigma_CX(E_rel) * v_rel
```

で CX rate を評価する。CX が実衝突として採択された場合、中性粒子速度は同じ sampled ion 速度 `vi` に置き換える。CX source も sampled `vi` を使い、運動量・エネルギー source は

```text
Sp_CX = dN * m_D * (v0 - vi)
Wi_CX = dN * (0.5 * m_D * |v0|^2 - 0.5 * m_D * |vi|^2)
```

で評価する。

`cx_model = 1` は、`src_read/monte/ntcros.f`, `ntfolw_atom.f`, `ntscor_cl.f`, `ntscor_tr.f` に寄せた比較用モデルである。CX rate は sampled ion 速度ではなく、背景 flow と ion 温度から作る有効相対速度を使う。

```text
v_rel^2 = |v0 - u_i|^2 + 8 * T_i * e / (pi * m_D)
E_rel   = 0.25 * m_D * v_rel^2 / e
R_CX    = n_i * sigma_CX_src_read(E_rel) * v_rel
```

CX 後の中性粒子速度は、`src_read` と同じく背景 flow 付き Maxwell 分布から新しく sample する。source は sampled `vi` ではなく、背景 flow と平均 thermal energy を使う。

```text
Sp_CX = dN * m_D * (v0 - u_i)
Wi_CX = dN * (0.5 * m_D * |v0|^2
              - 1.5 * T_i * e
              - 0.5 * m_D * |u_i|^2)
```

同じ seed と入力条件で `cx_model = 0` と `cx_model = 1` を切り替えると、現行 3D3V の sampled ion 実装と `src_read/monte` 寄せ実装の違いを比較できる。

## ビルドと実行

```bash
cd /home/mano/project/monte_simulations_ai/3d3v_ver1.0.0
make
make run
make plot
```

OpenMP で粒子ループを並列化する場合:

```bash
make OPENMP=1
OMP_NUM_THREADS=4 make OPENMP=1 run
```

計算速度を優先する場合は、bounds check なしの最適化ビルドも使える。

```bash
make OPENMP=1 FAST=1
OMP_NUM_THREADS=4 make OPENMP=1 FAST=1 run
```

逐次版と並列版は同じ粒子別 RNG を使うが、source の合算順は変わるため、浮動小数点の最終桁は完全一致しない場合がある。

主な出力は以下である。

- `run/sources.csv`: `CL`, `CL(avg)`, `TR` の `Sn/Sp/We/Wi`
- `run/source_stats.csv`: step 統計から計算した標準偏差・標準誤差・95% 信頼区間
- `run/balance.csv`: 粒子数、運動量、エネルギー収支
- `run/ntscrg.csv`: ion energy source の簡易 summary

`compare_estimator = .true.` の場合のみ、読み込み済み pretab table を `TR(pretab)` として追加出力する。

`make plot` は `run/delta_E_hist.csv` と `run/energy_hist.csv` から以下を作成する。

- `run/delta_E_hist.png`
- `run/energy_hist.png`

2つの実行結果を比較する場合は、標準では `run_01/delta_E_hist.csv` と
`run_02/delta_E_hist.csv` の `EL`, `CX`, `EL+CX` を1つのグラフに重ね描きする。

```bash
make plot-compare-deltaE
```

まだ `run_01` / `run_02` の配置にしていない場合は、入力ファイルを直接指定できる。

```bash
gnuplot -e "file1='run/delta_E_hist.csv'; file2='../other_run/delta_E_hist.csv'; label1='case A'; label2='case B'" plot_compare_deltaE.gp
```

## 変更履歴

1.0.0: 設計はDESIGN.md参照
1.0.1: gnuplotを追加、分散等も確認できるように
1.1.0: CXの比較を追加
