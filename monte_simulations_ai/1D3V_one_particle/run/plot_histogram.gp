# gnuplot script for collision histogram visualization
# Usage: gnuplot plot_histogram.gp

set terminal pngcairo enhanced size 1000,700 font 'Arial,12'
set output 'collision_histogram.png'

set title 'Energy Transfer Distribution: CX vs EL' font ',14'
set xlabel 'Energy Transfer ΔE [eV]' font ',12'
set ylabel 'Collision Frequency ν [s^{-1}]' font ',12'

set grid
set key right top font ',10'

# Log scale for y-axis
set logscale y

# Style settings
set style fill transparent solid 0.6 noborder

# Read metadata from file header
stats 'collision_histogram.dat' using 1 nooutput

# Plot both histograms
plot 'collision_histogram.dat' using 1:2 with boxes lc rgb '#3366CC' title 'Charge Exchange (CX)', \
     '' using 1:3 with boxes lc rgb '#DC3912' title 'Elastic Scattering (EL)', \
     0 with lines lc rgb 'black' lw 1.5 dt 2 title 'ΔE = 0'

print 'Plot saved to collision_histogram.png'

# Alternative: separate plots for CX and EL
set terminal pngcairo enhanced size 1200,500 font 'Arial,12'
set output 'collision_histogram_separate.png'

set multiplot layout 1,2

# CX histogram
set title 'Charge Exchange (CX)' font ',13'
set xlabel 'Energy Transfer ΔE [eV]'
set ylabel 'Reaction Rate [arb.]'
set logscale y
set grid
plot 'collision_histogram.dat' using 1:2 with boxes lc rgb '#3366CC' title 'CX rate', \
     0 with lines lc rgb 'black' lw 1.5 dt 2 notitle

# EL histogram
set title 'Elastic Scattering (EL)' font ',13'
plot 'collision_histogram.dat' using 1:3 with boxes lc rgb '#DC3912' title 'EL rate', \
     0 with lines lc rgb 'black' lw 1.5 dt 2 notitle

unset multiplot
print 'Plot saved to collision_histogram_separate.png'
