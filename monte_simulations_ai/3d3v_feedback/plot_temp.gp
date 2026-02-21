set terminal pngcairo size 800,600 enhanced font 'Verdana,10'
set output 'plasma_temp.png'

set datafile separator ','
set title 'Plasma Temperature and Density Evolution'
set xlabel 'Time [s]'

set multiplot layout 2,1

# First plot: Temperature
set ylabel 'Temperature [eV]'
#set logscale y
set format y "%.1e"
plot 'run/plasma_temp.csv' using 1:2 with lines lw 2 title 'T'

# Second plot: Density
set ylabel 'Electron Density [m^-3]'
#set logscale y
set format y "%.1e"
plot 'run/plasma_temp.csv' using 1:3 with lines lw 2 lc rgb 'red' title 'n_e'

unset multiplot
