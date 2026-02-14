#!/bin/bash
# E_neutral を 1〜10 eV でスキャンしてヒストグラムを生成

cd /home/mano/project/monte_simulations_ai/1D3V_one_particle/run

for E in 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0; do
    echo "Running E_neutral = ${E} eV..."
    
    # histogram_input.nml を更新
    cat > histogram_input.nml << EOF
&histogram_params
  T_neutral   = 1.0d0
  E_neutral   = ${E}d0
  T_ion       = 2.0d0
  n_bg        = 1.0d21
  n_samples   = 100000
  n_bins      = 80
  seed        = 54321
  cdf_file    = 'dd_00_elastic.cdf'
  output_file = 'collision_histogram.dat'
/
EOF
    
    # ヒストグラム計算を実行
    ./collision_histogram
    
    # 結果をhist/にコピー
    cp collision_histogram.dat ../hist/collision_histogram_E${E}eV.dat
done

echo "All done! Results saved to hist/"
