# 3D3V Monte Carlo Simulation ver 1.0.0

`project/monte_simulations_ai/3d3v_ver1.0.0/`はこれまでの非アナログモンテカルロシミュレーションをまとめたコードである。
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
`run_02/delta_E_hist.csv` を重ね描きする。

```bash
make plot-compare-deltaE
```

まだ `run_01` / `run_02` の配置にしていない場合は、入力ファイルを直接指定できる。

```bash
gnuplot -e "file1='run/delta_E_hist.csv'; file2='../other_run/delta_E_hist.csv'; label1='case A'; label2='case B'" plot_compare_deltaE.gp
```
