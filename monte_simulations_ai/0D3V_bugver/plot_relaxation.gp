#!/usr/bin/gnuplot
#===============================================================================
# plot_relaxation.gp
# エネルギー緩和過程の可視化
#===============================================================================

set terminal pngcairo size 800,600 enhanced font 'Arial,12'
set output 'relaxation.png'

set xlabel 'Time [s]'
set ylabel 'Energy [eV]'
set title '0D-3V Monte Carlo: Energy Relaxation'
set grid

# 目標エネルギー (3/2 * T_bg)
T_bg = 10.0
E_target = 1.5 * T_bg

# 目標線を描画
set arrow from graph 0,first E_target to graph 1,first E_target nohead lc rgb "red" dt 2 lw 2

# データプロット
plot 'run/relaxation.csv' using 1:2 with lines lw 2 lc rgb "blue" title 'Mean Energy', \
     E_target with lines lc rgb "red" dt 2 lw 2 title sprintf('Target (3/2 T_{bg} = %.1f eV)', E_target)
