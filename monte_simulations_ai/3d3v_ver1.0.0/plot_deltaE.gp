#!/usr/bin/env gnuplot
# Plot run/delta_E_hist.csv.

datafile = "run/delta_E_hist.csv"
outfile  = "run/delta_E_hist.png"

set datafile separator comma
set terminal pngcairo size 1200,800 enhanced font "Arial,12"
set output outfile

set title "Energy Transfer Histogram"
set xlabel "{/Symbol D}E [eV]"
set ylabel "Rate density [m^{-3} s^{-1}]"
set grid
set key outside right top
set xrange [-20:20]
set logscale y
set format y "10^{%L}"

set style line 1 lc rgb "#1f77b4" lw 2
set style line 2 lc rgb "#d62728" lw 2

plot datafile every ::1 using 1:($2 > 0 ? $2 : 1/0) with histeps ls 1 title "EL", \
     datafile every ::1 using 1:($3 > 0 ? $3 : 1/0) with histeps ls 2 title "CX"

print sprintf("Output: %s", outfile)
