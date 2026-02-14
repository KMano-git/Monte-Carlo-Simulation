# gnuplot script for particle trace visualization
# Usage: gnuplot plot_trace.gp

set terminal pngcairo enhanced size 1200,800 font 'Arial,12'
set output 'particle_trace.png'

set multiplot layout 2,2 title 'Single Particle Trace Analysis' font ',14'

# Load data
set datafile separator ','

# Panel 1: Energy vs Time
set xlabel 'Time [{/Symbol m}s]'
set ylabel 'Energy [eV]'
set title 'Particle Energy Evolution'
set grid
plot 'particle_trace.csv' using ($2*1e6):7 with steps lw 1.5 lc rgb 'blue' title 'E(t)', \
     '' using ($2*1e6):($8==1 ? $7 : 1/0) with points pt 7 ps 1.5 lc rgb 'red' title 'CX', \
     '' using ($2*1e6):($8==2 ? $7 : 1/0) with points pt 5 ps 1.5 lc rgb 'green' title 'EL', \
     '' using ($2*1e6):($8==3 ? $7 : 1/0) with points pt 4 ps 1.5 lc rgb 'purple' title 'EI'

# Panel 2: Position vs Time
set xlabel 'Time [{/Symbol m}s]'
set ylabel 'Position x [m]'
set title 'Particle Position'
set grid
plot 'particle_trace.csv' using ($2*1e6):3 with lines lw 1.5 lc rgb 'blue' title 'x(t)', \
     '' using ($2*1e6):($8==1 ? $3 : 1/0) with points pt 7 ps 1.5 lc rgb 'red' title 'CX', \
     '' using ($2*1e6):($8==2 ? $3 : 1/0) with points pt 5 ps 1.5 lc rgb 'green' title 'EL'

# Panel 3: Velocity components
set xlabel 'Time [{/Symbol m}s]'
set ylabel 'Velocity [m/s]'
set title 'Velocity Components'
set grid
plot 'particle_trace.csv' using ($2*1e6):4 with steps lw 1.5 lc rgb 'red' title 'v_x', \
     '' using ($2*1e6):5 with steps lw 1.5 lc rgb 'green' title 'v_y', \
     '' using ($2*1e6):6 with steps lw 1.5 lc rgb 'blue' title 'v_z'

# Panel 4: Energy change per collision (impulse-style bar chart)
set xlabel 'Time [{/Symbol m}s]'
set ylabel 'ΔE [eV]'
set title 'Energy Change at Each Collision'
set grid
set boxwidth 2e-2 absolute
set style fill solid 0.7

# Plot energy changes as impulses/bars, colored by collision type
# Filter only collision events (coll_type > 0)
plot 'particle_trace.csv' using ($2*1e6):($8==1 ? $9 : 1/0) with impulses lw 4 lc rgb 'red' title 'CX', \
     '' using ($2*1e6):($8==2 ? $9 : 1/0) with impulses lw 4 lc rgb 'forest-green' title 'EL', \
     '' using ($2*1e6):($8==3 ? $9 : 1/0) with impulses lw 4 lc rgb 'purple' title 'EI', \
     0 with lines lc rgb 'black' lw 0.5 dt 2 notitle

unset multiplot
print 'Plot saved to particle_trace.png'
