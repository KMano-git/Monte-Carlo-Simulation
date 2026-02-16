# Gnuplot script for 1D3V Flux Monte Carlo Simulation
# Usage: gnuplot plot_histogram.gp

set terminal pngcairo size 1200,800 enhanced font 'Arial,12'
set datafile separator ","

#==============================================================================
# Plot 1: 生存粒子数の時間変化
#==============================================================================
set output 'n_alive_vs_time.png'
set title '生存粒子数の時間変化' font ',14'
set xlabel 'Time [s]'
set ylabel 'N alive'
set grid
plot 'run/ntscrg.csv' using 1:2 with lines lw 2 title 'N alive' lc rgb '#0066CC'

#==============================================================================
# Plot 2: エネルギー交換の時間変化
#==============================================================================
set output 'energy_exchange_vs_time.png'
set title 'エネルギー交換量の時間変化' font ',14'
set xlabel 'Time [s]'
set ylabel 'ΔE [J]'
set grid
plot 'run/ntscrg.csv' using 1:3 with lines lw 2 title 'Total ΔE' lc rgb '#CC0000', \
     'run/ntscrg.csv' using 1:4 with lines lw 2 title 'CX ΔE' lc rgb '#00CC00', \
     'run/ntscrg.csv' using 1:5 with lines lw 2 title 'EL ΔE' lc rgb '#0000CC'

#==============================================================================
# Plot 3: エネルギーヒストグラム（最終ステップ）
#==============================================================================
set output 'energy_histogram.png'
set title 'エネルギー分布（各タイミング）' font ',14'
set xlabel 'Energy [eV]'
set ylabel 'Count'
set grid
set key top right

# エネルギーヒストグラムのCSV各行から読み取る
# 行番号でステップを区別
# 列: step, E_min, E_max, bin_width, n_alive, bin_1, bin_2, ...
# ビン中心 = E_min + (column_index - 5.5) * bin_width

# 個別のステップをプロットするにはawkで行を抽出
# 最後の行（最終ステップ）のエネルギーヒストグラム
stats 'run/energy_hist.csv' nooutput
N = STATS_records
set xrange [0:20]

plot for [row=2:N+1] 'run/energy_hist.csv' using (($2 + (column-5-0.5)*$4)):column every ::row-1::row-1 \
     with lines lw 2 title sprintf('row %d', row)

#==============================================================================
# Plot 4: 位置分布（最終ステップ）
#==============================================================================
set output 'position_distribution.png'
set title '位置分布（各タイミング）' font ',14'
set xlabel 'Position [m]'
set ylabel 'Count'
set grid
set autoscale x

stats 'run/statx.csv' nooutput
N = STATS_records

plot for [row=2:N+1] 'run/statx.csv' using (($2 + (column-5-0.5)*$4)):column every ::row-1::row-1 \
     with lines lw 2 title sprintf('row %d', row)

print "Plots generated successfully!"
