#!/usr/bin/env gnuplot
# Compare two delta_E_hist.csv files from 3d3v_ver1.0.1 runs.
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
if (!exists("outfile")) outfile = sprintf("delta_E_hist_compare_%s_%s.png", label1, label2)
if (!exists("xmin")) xmin = -20.0
if (!exists("xmax")) xmax = 20.0
if (!exists("ymin")) ymin = 1.0e15

set datafile separator comma
set terminal pngcairo size 1400,900 enhanced font "Arial,12"
set output outfile

set title "Delta E Histogram Comparison"
set grid
set key outside right top
set xrange [xmin:xmax]
set yrange [ymin:*]
set logscale y
set format y "10^{%L}"
set xlabel "{/Symbol D}E [eV]"
set ylabel "Rate [m^{-3} s^{-1}]"

# Use color for the quantity and line weight/brightness for the run.
# This stays readable with dense histogram steps better than dashed lines.
set style line 1 lc rgb "#005AB5" lw 4.0
set style line 2 lc rgb "#56B4E9" lw 2.0
set style line 3 lc rgb "#D55E00" lw 4.0
set style line 4 lc rgb "#E69F00" lw 2.0
set style line 5 lc rgb "#111111" lw 4.5
set style line 6 lc rgb "#7A7A7A" lw 2.2

plot '+' using (xmin):(ymin) with points pt 7 ps 0 notitle, \
     '+' using (xmin):(ymin*1.01) with points pt 7 ps 0 notitle, \
     file1 every ::1 using 1:($2 > 0 ? $2 : 1/0) with histeps ls 1 title sprintf("%s EL", label1), \
     file2 every ::1 using 1:($2 > 0 ? $2 : 1/0) with histeps ls 2 title sprintf("%s EL", label2), \
     file1 every ::1 using 1:($3 > 0 ? $3 : 1/0) with histeps ls 3 title sprintf("%s CX", label1), \
     file2 every ::1 using 1:($3 > 0 ? $3 : 1/0) with histeps ls 4 title sprintf("%s CX", label2), \
     file1 every ::1 using 1:(($2 + $3) > 0 ? ($2 + $3) : 1/0) with histeps ls 5 title sprintf("%s EL+CX", label1), \
     file2 every ::1 using 1:(($2 + $3) > 0 ? ($2 + $3) : 1/0) with histeps ls 6 title sprintf("%s EL+CX", label2)

unset output

print sprintf("Output: %s", outfile)
