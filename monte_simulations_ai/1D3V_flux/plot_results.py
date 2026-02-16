#!/usr/bin/env python3
"""1D3V Flux Monte Carlo Simulation — Results Plot"""
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['font.size'] = 12
plt.rcParams['figure.dpi'] = 150

# ========================================================================
# Read ntscrg.csv
# ========================================================================
time_arr, n_alive_arr = [], []
dE_total_arr, dE_cx_arr, dE_el_arr = [], [], []

with open('run/ntscrg.csv', 'r') as f:
    reader = csv.reader(f)
    next(reader)
    for row in reader:
        time_arr.append(float(row[0]))
        n_alive_arr.append(int(row[1]))
        dE_total_arr.append(float(row[2]))
        dE_cx_arr.append(float(row[3]))
        dE_el_arr.append(float(row[4]))

time_arr = np.array(time_arr)
n_alive_arr = np.array(n_alive_arr)
dE_total_arr = np.array(dE_total_arr)
dE_cx_arr = np.array(dE_cx_arr)
dE_el_arr = np.array(dE_el_arr)

# ========================================================================
# Read energy_hist.csv
# ========================================================================
energy_hist_data = []
with open('run/energy_hist.csv', 'r') as f:
    reader = csv.reader(f)
    next(reader)
    for row in reader:
        step = int(row[0])
        E_min = float(row[1])
        E_max = float(row[2])
        bin_width = float(row[3])
        n_alive = int(row[4])
        bins = [int(x) for x in row[5:]]
        energy_hist_data.append((step, E_min, E_max, bin_width, n_alive, bins))

# ========================================================================
# Read statx.csv
# ========================================================================
pos_hist_data = []
with open('run/statx.csv', 'r') as f:
    reader = csv.reader(f)
    next(reader)
    for row in reader:
        step = int(row[0])
        x_min = float(row[1])
        x_max = float(row[2])
        dx = float(row[3])
        n_alive = int(row[4])
        bins = [int(x) for x in row[5:]]
        pos_hist_data.append((step, x_min, x_max, dx, n_alive, bins))

# ========================================================================
# Plot
# ========================================================================
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle('1D3V Monte Carlo Simulation (Flux Model)', fontsize=16, fontweight='bold')

# --- (a) Alive particle count ---
ax = axes[0, 0]
ax.plot(time_arr * 1e6, n_alive_arr, color='#0066CC', linewidth=1.5)
ax.set_xlabel('Time [us]')
ax.set_ylabel('N alive')
ax.set_title('(a) Alive Particle Count vs Time')
ax.grid(True, alpha=0.3)

# --- (b) Energy exchange (50-step moving avg) ---
ax = axes[0, 1]
window = 50
if len(dE_total_arr) > window:
    dE_total_ma = np.convolve(dE_total_arr, np.ones(window)/window, mode='valid')
    dE_cx_ma = np.convolve(dE_cx_arr, np.ones(window)/window, mode='valid')
    dE_el_ma = np.convolve(dE_el_arr, np.ones(window)/window, mode='valid')
    t_ma = time_arr[window-1:]
    ax.plot(t_ma * 1e6, dE_total_ma, color='#CC0000', linewidth=1.5, label='Total')
    ax.plot(t_ma * 1e6, dE_cx_ma, color='#00AA00', linewidth=1.5, label='CX')
    ax.plot(t_ma * 1e6, dE_el_ma, color='#0000CC', linewidth=1.5, label='Elastic')
ax.set_xlabel('Time [us]')
ax.set_ylabel('Q [W/m3] (50-step avg)')
ax.set_title('(b) Energy Exchange vs Time')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

# --- (c) Energy distribution ---
ax = axes[1, 0]
colors = ['#AAAAAA', '#66AADD', '#3388BB', '#116699', '#CC3333', '#880000']
for idx, (step, E_min, E_max, bw, n_al, bins) in enumerate(energy_hist_data):
    if step == 0:
        continue
    n_bins = len(bins)
    E_centers = np.array([E_min + (i + 0.5) * bw for i in range(n_bins)])
    mask = E_centers < 20.0
    total = sum(bins)
    if total > 0:
        norm = np.array(bins, dtype=float) / total / bw
    else:
        norm = np.array(bins, dtype=float)
    c = colors[min(idx, len(colors)-1)]
    ax.plot(E_centers[mask], norm[mask], linewidth=1.5, color=c,
            label=f'step={step} (N={n_al})')
ax.set_xlabel('Energy [eV]')
ax.set_ylabel('Probability density [1/eV]')
ax.set_title('(c) Energy Distribution Evolution')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)
ax.set_xlim(0, 15)

# --- (d) Position distribution ---
ax = axes[1, 1]
for idx, (step, x_min, x_max, dx, n_al, bins) in enumerate(pos_hist_data):
    if step == 0:
        continue
    n_bins = len(bins)
    x_centers = np.array([x_min + (i + 0.5) * dx for i in range(n_bins)])
    c = colors[min(idx, len(colors)-1)]
    ax.plot(x_centers * 100, bins, linewidth=1.5, color=c,
            label=f'step={step} (N={n_al})')
ax.set_xlabel('Position [cm]')
ax.set_ylabel('Count')
ax.set_title('(d) Position Distribution Evolution')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('simulation_results.png', bbox_inches='tight')
print("Saved: simulation_results.png")
