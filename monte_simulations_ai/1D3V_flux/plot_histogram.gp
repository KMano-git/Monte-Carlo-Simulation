# Gnuplot script for 1D3V Flux Monte Carlo Simulation
# Usage:
#   bash preprocess.sh      (CSVを2列形式に変換)
#   gnuplot plot_histogram.gp
#
# 4パネル: (a)粒子数, (b)エネルギー交換, (c)エネルギー分布, (d)位置分布

set terminal pngcairo size 1400,1000 enhanced font 'Arial,11'
set output 'simulation_results.png'
set datafile separator ","

set multiplot layout 2,2 title '1D3V Monte Carlo Simulation (Flux Model)' font ',16'

#==============================================================================
# (a) Alive Particle Count vs Time
#==============================================================================
set title '(a) Alive Particle Count vs Time' font ',13'
set xlabel 'Time [us]'
set ylabel 'N alive'
set grid
set key off
plot 'run/ntscrg.csv' using ($1*1e6):2 every ::1 with lines lw 2 lc rgb '#0066CC'

#==============================================================================
# (b) Energy Exchange vs Time
#==============================================================================
set title '(b) Energy Exchange vs Time' font ',13'
set xlabel 'Time [us]'
set ylabel 'dE [J]'
set grid
set key top right
plot 'run/ntscrg.csv' using ($1*1e6):3 every ::1 smooth csplines lw 2 lc rgb '#CC0000' title 'Total', \
     'run/ntscrg.csv' using ($1*1e6):4 every ::1 smooth csplines lw 2 lc rgb '#00AA00' title 'CX', \
     'run/ntscrg.csv' using ($1*1e6):5 every ::1 smooth csplines lw 2 lc rgb '#0000CC' title 'Elastic'

#==============================================================================
# (c) Energy Distribution Evolution
# preprocess.sh で生成した /tmp/ehist_*.dat を使用
#==============================================================================
set title '(c) Energy Distribution Evolution' font ',13'
set xlabel 'Energy [eV]'
set ylabel 'Count'
set grid
set key top right
set xrange [0:15]
set autoscale y
set datafile separator whitespace

plot '/tmp/ehist_2.dat' using 1:2 with lines lw 2 lc rgb '#66AADD' title 'step=100', \
     '/tmp/ehist_3.dat' using 1:2 with lines lw 2 lc rgb '#3388BB' title 'step=500', \
     '/tmp/ehist_4.dat' using 1:2 with lines lw 2 lc rgb '#116699' title 'step=1000', \
     '/tmp/ehist_5.dat' using 1:2 with lines lw 2 lc rgb '#CC3333' title 'step=2000', \
     '/tmp/ehist_6.dat' using 1:2 with lines lw 2 lc rgb '#880000' title 'step=5000'

#==============================================================================
# (d) Position Distribution Evolution
# preprocess.sh で生成した /tmp/phist_*.dat を使用
#==============================================================================
set title '(d) Position Distribution Evolution' font ',13'
set xlabel 'Position [cm]'
set ylabel 'Count'
set grid
set key top right
set autoscale x
set autoscale y

plot '/tmp/phist_2.dat' using 1:2 with lines lw 2 lc rgb '#66AADD' title 'step=100', \
     '/tmp/phist_3.dat' using 1:2 with lines lw 2 lc rgb '#3388BB' title 'step=500', \
     '/tmp/phist_4.dat' using 1:2 with lines lw 2 lc rgb '#116699' title 'step=1000', \
     '/tmp/phist_5.dat' using 1:2 with lines lw 2 lc rgb '#CC3333' title 'step=2000', \
     '/tmp/phist_6.dat' using 1:2 with lines lw 2 lc rgb '#880000' title 'step=5000'

unset multiplot
print "Saved: simulation_results.png"
