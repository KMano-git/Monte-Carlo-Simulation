#!/usr/bin/gnuplot
#===============================================================================
# plot_histogram.gp
# energy_hist.csv からエネルギー分布ヒストグラムをプロット
# usage: cd 0D3V/ && gnuplot plot_histogram.gp
#===============================================================================

set terminal pngcairo size 1400,900 enhanced font 'Arial,12'
set output 'energy_histogram.png'

# パラメータ (input.nml と合わせる)
E_min_val  = 0.0
E_max_val  = 150.0
n_bins     = 400
bin_width  = (E_max_val - E_min_val) / n_bins
T_i        = 2.0   # 背景イオン温度 [eV]

# Maxwell分布 (3D速度空間)
maxwell(E) = 2.0 * sqrt(E/pi) / (T_i**1.5) * exp(-E/T_i)

# 前処理: 各行を個別ファイルに展開
system("awk -F',' 'NR>1 && NF>5 { \
    total=0; for(i=6;i<=NF;i++) total+=$i; \
    if(total==0) next; \
    bw=$4+0; emin=$2+0; \
    fname=\"/tmp/hist_row_\" (NR-2) \".dat\"; \
    system(\"rm -f \" fname); \
    for(i=6;i<=NF;i++) print emin+(i-6+0.5)*bw, ($i+0)/(total*bw) > fname; \
    close(fname) \
}' run/energy_hist.csv")

# データ行数を取得
n_rows = int(system("awk -F',' 'NR>1 && NF>5{n++} END{print n}' run/energy_hist.csv"))

# 各行のステップ番号を配列に格納
do for [row=0:n_rows-1] {
    eval sprintf("step_%d = '"."%s'", row, \
        system(sprintf("awk -F',' 'NR==%d{gsub(/ /,\"\",$1); printf $1}' run/energy_hist.csv", row+2)))
}

# 表示設定
set xrange [0:20]
set xlabel 'Energy [eV]'
set ylabel 'Probability Density [1/eV]'
set grid
set style fill solid 0.3
set boxwidth bin_width

set multiplot layout 2,3 title '0D-3V Monte Carlo: Energy Distribution Evolution' font ',14'

do for [row=0:n_rows-1] {
    eval sprintf("stitle = step_%d", row)
    set title sprintf('Step %s (dt=10ns)', stitle)

    tmpfile = sprintf('/tmp/hist_row_%d.dat', row)

    plot tmpfile using 1:2 with boxes lc rgb "blue" notitle, \
         maxwell(x) with lines lw 2 lc rgb "red" title sprintf('Maxwell (T=%.1f eV)', T_i)
}

unset multiplot

# 一時ファイル削除
system("rm -f /tmp/hist_row_*.dat")

print "Output: energy_histogram.png"
