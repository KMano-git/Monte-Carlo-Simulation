#!/usr/bin/gnuplot
#===============================================================================
# plot_histogram.gp
# エネルギー分布とMaxwell分布の比較
#===============================================================================

set terminal pngcairo size 800,600 enhanced font 'Arial,12'
set output 'energy_distribution.png'

set xlabel 'Energy [eV]'
set ylabel 'Probability Density'
set title '0D-3V Monte Carlo: Final Energy Distribution vs Maxwell'
set grid

# Maxwell分布の理論式 (3D)
T_bg = 10.0
maxwell(E) = 2.0*sqrt(E/pi)/(T_bg**1.5) * exp(-E/T_bg)

# データプロット
plot 'run/energy_dist.csv' using 1:2 with boxes fill solid 0.3 lc rgb "blue" title 'Simulation', \
     maxwell(x) with lines lw 2 lc rgb "red" title sprintf('Maxwell (T=%.1f eV)', T_bg)
