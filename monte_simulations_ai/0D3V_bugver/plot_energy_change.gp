set terminal pngcairo size 800,600 enhanced font 'Arial,12'
set output 'energy_change.png'

set xlabel 'Time [s]'
set ylabel 'Average Energy Change per Particle [eV]'
set title 'Energy Change by Collision Type'
set grid

# 凡例の位置
set key top right

# 0のライン
set xzeroaxis lt 1 lc rgb "black"

plot 'run/ntscrg.csv' using 1:2 with lines lw 2 lc rgb "black" title 'Total', \
     'run/ntscrg.csv' using 1:3 with lines lw 2 lc rgb "blue" title 'Charge Exchange', \
     'run/ntscrg.csv' using 1:4 with lines lw 2 lc rgb "red" title 'Elastic Scattering'
