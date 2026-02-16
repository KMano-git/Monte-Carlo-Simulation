#!/bin/bash
# energy_hist.csv / statx.csv を gnuplot 用の2列形式に変換する
# Usage: bash preprocess.sh
# 出力: /tmp/ehist_*.dat, /tmp/phist_*.dat

# エネルギーヒストグラム
row=0
tail -n +2 run/energy_hist.csv | while IFS=, read -r line; do
    row=$((row+1))
    step=$(echo "$line" | cut -d, -f1)
    emin=$(echo "$line" | cut -d, -f2)
    bw=$(echo "$line" | cut -d, -f4)
    nalive=$(echo "$line" | cut -d, -f5)
    # step=0はスキップ
    if [ "$step" = "0" ]; then continue; fi
    outfile="/tmp/ehist_${row}.dat"
    echo "# step=$step N=$nalive" > "$outfile"
    # ビン値を1列ずつ出力
    nbins=$(echo "$line" | awk -F, '{print NF-5}')
    for i in $(seq 1 $nbins); do
        col=$((i+5))
        center=$(echo "$emin $bw $i" | awk '{printf "%.4f", $1+($3-0.5)*$2}')
        val=$(echo "$line" | cut -d, -f$col)
        echo "$center $val" >> "$outfile"
    done
done

# 位置ヒストグラム
row=0
tail -n +2 run/statx.csv | while IFS=, read -r line; do
    row=$((row+1))
    step=$(echo "$line" | cut -d, -f1)
    xmin=$(echo "$line" | cut -d, -f2)
    dx=$(echo "$line" | cut -d, -f4)
    nalive=$(echo "$line" | cut -d, -f5)
    if [ "$step" = "0" ]; then continue; fi
    outfile="/tmp/phist_${row}.dat"
    echo "# step=$step N=$nalive" > "$outfile"
    nbins=$(echo "$line" | awk -F, '{print NF-5}')
    for i in $(seq 1 $nbins); do
        col=$((i+5))
        center=$(echo "$xmin $dx $i" | awk '{printf "%.6f", ($1+($3-0.5)*$2)*100}')
        val=$(echo "$line" | cut -d, -f$col)
        echo "$center $val" >> "$outfile"
    done
done

echo "Preprocessed data files written to /tmp/ehist_*.dat and /tmp/phist_*.dat"
