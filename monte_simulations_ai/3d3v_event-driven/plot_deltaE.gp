# plot_deltaE.gp
# Plot the energy transfer (delta E) histograms for EL and CX from deltaE_hist.csv

# Output settings
set terminal pngcairo size 1000,800 font "Arial,12"
set output 'deltaE_histogram.png'

# General plot settings
set title "Energy Transfer (Delta E) Histogram" font ",14"
set xlabel "Delta E [eV]" font ",12"
set ylabel "Reaction Rate Density [events / (m^3 * s)]" font ",12"
set grid

# Log scale for y-axis to see both large and small variations
set logscale y
set format y "10^{%L}"
set yrange [1e15:*] # 最小値を調整して見やすくする

# X-axis range (adjust if necessary, but leaving it auto usually shows the full -100 to 100 eV range)
set xrange [-20:20]

# Set line styles
set style line 1 lc rgb '#0060ad' pt 7 ps 0.5 lt 1 lw 2 # Blue for EL
set style line 2 lc rgb '#dd181f' pt 7 ps 0.5 lt 1 lw 2 # Red for CX

# We plot the last block in the CSV file, which corresponds to the final output step
# The data starts at line 3 of the final block
plot 'run/deltaE_hist.csv' index (system("grep -c '^step' run/deltaE_hist.csv")-1) every ::1 using 1:2 with histeps ls 1 title 'Elastic (EL)', \
     'run/deltaE_hist.csv' index (system("grep -c '^step' run/deltaE_hist.csv")-1) every ::1 using 1:3 with histeps ls 2 title 'Charge Exchange (CX)'

print "Generated deltaE_histogram.png"
