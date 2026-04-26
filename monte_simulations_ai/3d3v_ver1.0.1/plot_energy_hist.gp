#!/usr/bin/env gnuplot
# Plot run/energy_hist.csv. The CSV stores one timestep per row and one
# histogram bin per column, so this script reshapes it with awk first.

datafile = "run/energy_hist.csv"
longfile = "run/energy_hist_long.dat"
outfile  = "run/energy_hist.png"

# Keep these consistent with &diagnostics in run/input.nml.
E_min = 0.0
E_max = 150.0
n_bins = 400
bin_width = (E_max - E_min) / n_bins

system(sprintf("awk -F',' -v emin=%g -v bw=%g 'NR > 1 { step=$1+0; total=$3+0; for (i=4; i<=NF; i++) { e=emin+(i-4+0.5)*bw; y=(total>0 ? ($i+0)/(total*bw) : 0); print step, e, y } print \"\"; print \"\" }' %s > %s", E_min, bin_width, datafile, longfile))

n_steps = int(system(sprintf("awk -F',' 'NR > 1 {n++} END{print n+0}' %s", datafile)))

set terminal pngcairo size 1400,900 enhanced font "Arial,12"
set output outfile

set title "Neutral Energy Distribution"
set xlabel "Energy [eV]"
set ylabel "Probability density [eV^{-1}]"
set grid
set key outside right top
set xrange [0:30]
set yrange [0:*]

set style line 1 lc rgb "#1f77b4" lw 2
set style line 2 lc rgb "#ff7f0e" lw 2
set style line 3 lc rgb "#2ca02c" lw 2
set style line 4 lc rgb "#d62728" lw 2
set style line 5 lc rgb "#9467bd" lw 2
set style line 6 lc rgb "#8c564b" lw 2

plot for [i=0:n_steps-1] longfile index i using 2:3 with lines lw 2 title sprintf("step %s", system(sprintf("awk -F',' 'NR==%d {print $1; exit}' %s", i+2, datafile)))

print sprintf("Output: %s", outfile)
system(sprintf("rm -f %s", longfile))
