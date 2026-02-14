# gnuplot script for plotting CX-EL difference
# Usage: gnuplot plot_cx_el_diff.gp

set terminal pngcairo enhanced size 1000,1200 font 'Arial,11'
set output 'hist/cx_el_diff.png'

set multiplot layout 2,1 title 'CX - EL Collision Frequency Difference (T_{ion} = 2.0 eV)' font ',14'

# 10色のカラーパレット（区別しやすい色）
set linetype 1 lc rgb '#1b0c9e' # 濃い青
set linetype 2 lc rgb '#0072B2' # 青
set linetype 3 lc rgb '#009E73' # 緑
set linetype 4 lc rgb '#56B4E9' # 水色
set linetype 5 lc rgb '#F0E442' # 黄
set linetype 6 lc rgb '#E69F00' # オレンジ
set linetype 7 lc rgb '#D55E00' # 赤オレンジ
set linetype 8 lc rgb '#CC0033' # 赤
set linetype 9 lc rgb '#882288' # 紫
set linetype 10 lc rgb '#AA4499' # ピンク紫

# Common settings
set xlabel 'ΔE [eV]'
set ylabel 'Frequency Difference [s^{-1}/eV]'
set grid
set key outside right top
set xrange [-12:8]

# Panel 1: Positive region (where CX > EL)
set title 'CX - EL Difference'
set yrange [-1e6:2e6]
plot 'hist/collision_histogram_E1.0eV.dat' using 1:($2-$3) with lines lw 2 lt 1 title 'E_n = 1 eV', \
     'hist/collision_histogram_E2.0eV.dat' using 1:($2-$3) with lines lw 2 lt 2 title 'E_n = 2 eV', \
     'hist/collision_histogram_E3.0eV.dat' using 1:($2-$3) with lines lw 2 lt 3 title 'E_n = 3 eV', \
     'hist/collision_histogram_E4.0eV.dat' using 1:($2-$3) with lines lw 2 lt 4 title 'E_n = 4 eV', \
     'hist/collision_histogram_E5.0eV.dat' using 1:($2-$3) with lines lw 2 lt 5 title 'E_n = 5 eV', \
     'hist/collision_histogram_E6.0eV.dat' using 1:($2-$3) with lines lw 2 lt 6 title 'E_n = 6 eV', \
     'hist/collision_histogram_E7.0eV.dat' using 1:($2-$3) with lines lw 2 lt 7 title 'E_n = 7 eV', \
     'hist/collision_histogram_E8.0eV.dat' using 1:($2-$3) with lines lw 2 lt 8 title 'E_n = 8 eV', \
     'hist/collision_histogram_E9.0eV.dat' using 1:($2-$3) with lines lw 2 lt 9 title 'E_n = 9 eV', \
     'hist/collision_histogram_E10.0eV.dat' using 1:($2-$3) with lines lw 2 lt 10 title 'E_n = 10 eV'

# Panel 2: Log scale absolute difference  
set title '|CX - EL| Absolute Difference (Log Scale)'
set ylabel '|Frequency Difference| [s^{-1}/eV]'
set logscale y
set yrange [1e3:1e8]
plot 'hist/collision_histogram_E1.0eV.dat' using 1:(abs($2-$3)) with lines lw 2 lt 1 title 'E_n = 1 eV', \
     'hist/collision_histogram_E2.0eV.dat' using 1:(abs($2-$3)) with lines lw 2 lt 2 title 'E_n = 2 eV', \
     'hist/collision_histogram_E3.0eV.dat' using 1:(abs($2-$3)) with lines lw 2 lt 3 title 'E_n = 3 eV', \
     'hist/collision_histogram_E4.0eV.dat' using 1:(abs($2-$3)) with lines lw 2 lt 4 title 'E_n = 4 eV', \
     'hist/collision_histogram_E5.0eV.dat' using 1:(abs($2-$3)) with lines lw 2 lt 5 title 'E_n = 5 eV', \
     'hist/collision_histogram_E6.0eV.dat' using 1:(abs($2-$3)) with lines lw 2 lt 6 title 'E_n = 6 eV', \
     'hist/collision_histogram_E7.0eV.dat' using 1:(abs($2-$3)) with lines lw 2 lt 7 title 'E_n = 7 eV', \
     'hist/collision_histogram_E8.0eV.dat' using 1:(abs($2-$3)) with lines lw 2 lt 8 title 'E_n = 8 eV', \
     'hist/collision_histogram_E9.0eV.dat' using 1:(abs($2-$3)) with lines lw 2 lt 9 title 'E_n = 9 eV', \
     'hist/collision_histogram_E10.0eV.dat' using 1:(abs($2-$3)) with lines lw 2 lt 10 title 'E_n = 10 eV'

unset multiplot
print 'Plot saved to hist/cx_el_diff.png'
