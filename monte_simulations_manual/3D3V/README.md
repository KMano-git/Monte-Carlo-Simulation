# 3D3V Monte Carlo Simulation (実質0次元)

核融合周辺プラズマにおける中性粒子輸送・挙動解析に用いられるNEUT2Dの機能を段階的に実装し、検証するためのコード。
非アナログ法、イベント駆動、collision/track-length estimatorで実装されている。

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


