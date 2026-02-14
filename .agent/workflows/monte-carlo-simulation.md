---
description: Fortranを使用したシミュレーションコード開発。具体的な解析対象は核融合周辺プラズマであり、特に中性粒子をモンテカルロ法で詳しく追う。
---

# Role
あなたはFortranによる計算物理シミュレーションの専門家です。
数値計算に強く、モダンなFortran (Fortran 2008 standard) のコーディング規約を遵守します。

# Goal
背景プラズマが固定された環境下での、中性粒子の振る舞いを解析するための「1次元 モンテカルロシミュレーション (Test Particle Monte Carlo)」を作成してください。

# 物理モデルの仕様
1. **次元**: 空間1次元 ($x$)、速度3次元 ($v_x, v_y, v_z$)
2. **時間発展**: 一定の時間刻み $\Delta t$ で粒子を推進させる (Time-driven)。
3. **衝突判定アルゴリズム (No Null-Collision)**:
    - ステップごとに全衝突頻度 $\nu_{total} = v \cdot n_{bg} \cdot (\sigma_{CX} + \sigma_{EL} + \sigma_{EI})$ を計算。
    - 衝突確率 $P = 1 - \exp(-\nu_{total} \Delta t)$ と乱数 $R$ を比較して衝突有無を決定。
    - 衝突する場合、断面積の比率に応じてプロセス（CX か Elastic か EI）を選択。
4. **対象プロセス**:
    - **荷電交換 (Charge Exchange, CX)**: 粒子が背景粒子とエネルギー・運動量を交換する（速度が背景温度の分布に入れ替わる）。
    - **弾性衝突 (Elastic Scattering)**: エネルギー保存、運動量の方向散乱。衝突断面積と散乱角データは"dd_00_elastic.cdf"から読み込む。
    - **電離衝突**: 粒子が背景電子と衝突し、イオン化する反応。閾値エネルギー以上のエネルギーを受け取った場合にのみ発生する。二次電子は追跡しません。
5. **背景**: 固定された密度一定のプラズマ ($n_{bg}$)、温度 $T_{bg}$。
6. **外場**: 一定の電磁場を実装（今のところは電場も磁場も0で考える）。

# 重要: スコアリング手法の実装 (Scoring Methods)
2つの手法を実装し、モジュール `scoring` 内で管理してください。両方の結果を同時に計算し、比較できるようにします。

## 1. Collision Estimator (CL) - `score_cl`
- **タイミング**: 乱数判定で「衝突が発生した (Collided)」時のみ呼び出す。
- **ロジック**:
    - 選択されたプロセス ($k$) に応じたエネルギー損失/変化量 $\Delta E$ を記録。
    - 参考: `Grid_Energy(ix) += Weight * Delta_E`

## 2. Track-Length Estimator (TL) - `score_tr`
- **タイミング**: 衝突の有無に関わらず、毎ステップ ($\Delta t$) 呼び出す。
- **ロジック**:
    - そのステップでの期待値（確率的な寄与）を積分する。
    - `Contribution = Weight * (Rate_CX * Loss_CX + Rate_EL * Loss_EL) * dt`
    - ここで `Rate = n_{bg} * sigma * v`。

# コーディング要件 (Modern Fortran)
- **ファイル構成**: 単一ファイル (`main.f90`) でも可だが、内部で `module` を使って整理すること。
- **モジュール構成**:
    - `module types`: 粒子 (`type particle`) やグリッドの定義。
    - `module cross_sections`: $\sigma_{CX}(E)$ と $\sigma_{EL}(E)$ のモデル関数（定数または $1/\sqrt{E}$ 簡易モデル）。
    - `module scoring`: `score_cl`, `score_tr` サブルーチン。
    - `module dynamics`: 粒子の推進と衝突判定分岐。
- **出力**:
    - `results.csv`: 位置 $x$ ごとの「CL法によるエネルギー付与」と「TL法によるエネルギー付与」を列挙して比較可能にする。

# Workflow Steps

## Step 1: 設計 (Design)
- 衝突判定のロジック（分岐比の計算）を数式で示してください。
- `ntscor_cl.f` (実衝突時の処理) と `ntscor_tr.f` (飛行中の積分) の概念を、今回の単純化モデルにどう適用するか整理してください。
- 外部ファイル読み込みとデータ補間のロジックを設計してください。1次補完で十分です。

## Step 2: 実装 (Implementation)
- 設計に基づきFortranコードを記述してください。
- 変数宣言は必ず `implicit none` を使用してください。
- 定数は `double precision` (kind=8) を使用してください。
- メインループでは、粒子が領域外にでた場合の境界条件も考慮してください。 

## Step 3: コンパイル・実行ガイド
- `gfortran` でのコンパイルコマンドと、出力データのgnuplotでの描画例を示してください。