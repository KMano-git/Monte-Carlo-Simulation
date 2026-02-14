# モンテカルロ中性粒子輸送シミュレーション

1次元空間・3次元速度空間における中性粒子（重水素）の輸送をモンテカルロ法でシミュレートするFortranコードです。

## 特徴
- **物理過程**: 荷電交換(CX)、弾性散乱(EL)、電子衝突電離(EI)
- **スコアリング**: Collision Estimator (CL法) と Track-Length Estimator (TL法) の併用
- **出力**: 物理単位（W/m², m⁻³）でのパワー密度および粒子密度分布

## ディレクトリ構成
```
project/
├── code/              # ソースコード
│   ├── constants.f90      # 物理定数
│   ├── types.f90          # データ型定義
│   ├── cdf_reader.f90     # CDFファイル読み込み
│   ├── cross_sections.f90 # 断面積計算
│   ├── scoring.f90        # スコアリング
│   ├── dynamics.f90       # 粒子推進・衝突
│   ├── io.f90             # 入出力処理
│   └── main.f90           # メインプログラム
├── run/               # 実行ディレクトリ
│   ├── input.nml          # 入力パラメータファイル
│   └── dd_00_elastic.cdf  # 断面積データファイル
├── Makefile           # ビルドファイル
└── README.md          # 本マニュアル
```

## コンパイルと実行
```bash
make
cd run
./monte_carlo
```

## 入力パラメータ (input.nml)
```fortran
&simulation
  n_particles = 10000    ! テスト粒子数
  n_steps     = 10000    ! ステップ数
  n_grid      = 100      ! グリッド数
  dt          = 1.0d-8   ! 時間刻み [s]
  x_min       = 0.0d0    ! 領域下限 [m]
  x_max       = 0.1d0    ! 領域上限 [m]
  output_file = 'ntscrg.csv' ! パワー出力ファイル
/

&plasma_nml
  n_bg = 1.0d19  ! 背景イオン密度 [m^-3]
  T_bg = 2.0d0   ! 背景イオン温度 [eV]
  T_e  = 2.0d0   ! 電子温度 [eV]
/

&particle_init
  gamma_in = 1.0d23  ! 入射粒子フラックス [m^-2 s^-1]
  E_init   = 5.0d0   ! 初期エネルギー [eV]
/
```

## 出力ファイル
- **ntscrg.csv**: 各反応（CX, EL, EI）によるパワー密度分布 [W/m²]
- **ntden.csv**: 中性粒子密度分布 [m⁻³]（フラックス密度換算）
