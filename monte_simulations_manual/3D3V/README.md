# 3D3V Monte Carlo Simulation (実質0次元)

核融合周辺プラズマにおける中性粒子輸送・挙動解析に用いられるNEUT2Dの機能を段階的に実装し、検証するためのコード。
非アナログ法、イベント駆動、collision/track-length estimatorで実装されている。
AI生成されたコードをもとに、色々と手直しを含めて作成する最終的な完成コードにする予定。

## 参考にするコードについて

 - ntflow_atom.f: main.f90に一番近い形。EIの重み減衰、CX/ELの選択、scoringのタイミング、衝突後の更新について
 - ntscor_cl.f: Collision estimator のCXとEIの参考に
 - ntscor_tr.f: Track-length estimator について飛跡長ごとにどのような処理、平均化をしているか
 - ntcros.f: 反応率 srct をどう計算しているかの本体。荷電交換の扱いを見るならここ
 - ntvel_el.f: ELの実衝突サンプリングと、CL用の pcl_eng, pcl_mvp などを計算する
 - ntsctr_el.f: Track-Length における弾性衝突の積分カーネルを組み立てる本体
 - celcom.f: 共通変数定義
 - ntflow.f: 原子分子共通の処理。TRを呼ぶタイミング
 - ntmont.f: 上位ドライバ

## 満たすべき条件

衝突判定: Null collision 法
時間ステップ: Event-driven
スコアリング方法: Analog/ Collision/ Track-length
    平均化: テーブル参照法


## ディレクトリ構成

## 実行方法
1. `make` or `make OPENMP=1` でコンパイル
2. `input.nml` でパラメータを設定
3. `./monte_carlo*` で実行

## 使用している断面積データについて

現在は `dd_00_elastic.cdf` (Bachmann(1995)) のデータを元に計算している

### `dd_00_elastic.cdf` のデータ構造

`dd_00_elastic.cdf` は、以下の順序で配列データが格納されている（すべて連続した1次元配列として並んでいる）。

1. `cross_section` (101個)
2. `reaction_rate` (2601個)
3. `I_1_0` (2601個)
4. `I_1_1*up` (2601個)
5. `I_1_2*up^2` (2601個)
6. `scattering_angle` (12801個)
7. `sigv_max` (1個)
8. `angle_min` (1個)


