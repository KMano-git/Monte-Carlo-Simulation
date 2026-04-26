#!/usr/bin/env gnuplot
# Compare two delta_E_hist.csv files from 3d3v_ver1.0.0 runs.
#
# Default:
#   gnuplot plot_compare_deltaE.gp
#
# Override paths/labels while the run_01/run_02 directory layout is not ready:
#   gnuplot -e "file1='run/delta_E_hist.csv'; file2='../other_run/delta_E_hist.csv'; label1='case A'; label2='case B'" plot_compare_deltaE.gp

if (!exists("file1")) file1 = "run_01/delta_E_hist.csv"
if (!exists("file2")) file2 = "run_02/delta_E_hist.csv"
if (!exists("label1")) label1 = "run_01"
if (!exists("label2")) label2 = "run_02"
if (!exists("outfile")) outfile = "delta_E_hist_compare.png"
if (!exists("xmin")) xmin = -20.0
if (!exists("xmax")) xmax = 20.0
if (!exists("ymin")) ymin = 1.0e15

set datafile separator comma
set terminal pngcairo size 1400,1000 enhanced font "Arial,12"
set output outfile

set grid
set key outside right top
set xrange [xmin:xmax]
set yrange [ymin:*]
set logscale y
set format y "10^{%L}"
set xlabel "{/Symbol D}E [eV]"

set style line 1 lc rgb "#1f77b4" lw 2
set style line 2 lc rgb "#ff7f0e" lw 2 dt 2
set style line 3 lc rgb "#2ca02c" lw 2
set style line 4 lc rgb "#d62728" lw 2 dt 2

set multiplot layout 2,1 title "Delta E Histogram Comparison"

set title "Elastic (EL)"
set ylabel "EL rate [m^{-3} s^{-1}]"
plot '+' using (xmin):(ymin) with points pt 7 ps 0 notitle, \
     '+' using (xmin):(ymin*1.01) with points pt 7 ps 0 notitle, \
     file1 every ::1 using 1:($2 > 0 ? $2 : 1/0) with histeps ls 1 title sprintf("%s EL", label1), \
     file2 every ::1 using 1:($2 > 0 ? $2 : 1/0) with histeps ls 2 title sprintf("%s EL", label2)

set title "Charge Exchange (CX)"
set ylabel "CX rate [m^{-3} s^{-1}]"
plot '+' using (xmin):(ymin) with points pt 7 ps 0 notitle, \
     '+' using (xmin):(ymin*1.01) with points pt 7 ps 0 notitle, \
     file1 every ::1 using 1:($3 > 0 ? $3 : 1/0) with histeps ls 3 title sprintf("%s CX", label1), \
     file2 every ::1 using 1:($3 > 0 ? $3 : 1/0) with histeps ls 4 title sprintf("%s CX", label2)

unset multiplot
unset output

print sprintf("Output: %s", outfile)
