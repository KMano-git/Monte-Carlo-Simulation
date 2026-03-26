import os
import re
import sys
import numpy as np
import scipy.integrate as integrate
from scipy.interpolate import interp1d

print("Starting calc_I_kernel.py...", flush=True)

# File paths
DEFAULT_IN_CDF = '../test/cdf_change/dd_00_elastic_pure_el_angle.cdf'
DEFAULT_OUT_CDF = 'dd_00_elastic_pure_el_angle_fixed.cdf'

if len(sys.argv) >= 2:
    IN_CDF = sys.argv[1]
else:
    IN_CDF = DEFAULT_IN_CDF

if len(sys.argv) >= 3:
    OUT_CDF = sys.argv[2]
else:
    OUT_CDF = DEFAULT_OUT_CDF

def parse_array(text, name):
    pat = re.compile(rf'\b{re.escape(name)}\s*=\s*(.*?)\s*;', re.DOTALL)
    m = pat.search(text)
    if not m: return []
    raw = m.group(1).replace('\n', ' ')
    vars = [x.strip() for x in raw.split(',') if x.strip()]
    return vars

def parse_int_array(text, name):
    return [int(v) for v in parse_array(text, name)]
def parse_float_array(text, name):
    return [float(v) for v in parse_array(text, name)]
def make_axis(n, vmin, vmax, spacing):
    spacing = spacing.replace('"', '').strip()
    if spacing == 'log' and vmin > 0 and vmax > 0:
        return np.logspace(np.log10(vmin), np.log10(vmax), n)
    else:
        return np.linspace(vmin, vmax, n)

print("Reading input CDF...", flush=True)
with open(IN_CDF, 'r') as f:
    orig_text = f.read()

xs_data_tab_str = parse_array(orig_text, 'xs_data_tab')
xs_data_tab = np.array([float(x) for x in xs_data_tab_str])

xs_data_base = parse_int_array(orig_text, 'xs_data_base')
xs_data_inc  = parse_int_array(orig_text, 'xs_data_inc')

min_f = parse_float_array(orig_text, 'xs_min'); xs_min = [min_f[i*3:(i+1)*3] for i in range(20)]
spc_r = parse_array(orig_text, 'xs_spacing'); xs_spacing = [spc_r[i*4:(i+1)*4] for i in range(20)]
tab_f = parse_int_array(orig_text, 'xs_tab_index'); xs_tab_index = [tab_f[i*3:(i+1)*3] for i in range(20)]
max_f = parse_float_array(orig_text, 'xs_max'); xs_max = [max_f[i*3:(i+1)*3] for i in range(20)]

def get_data(idx):
    return np.array(xs_data_tab[xs_data_base[idx] : xs_data_base[idx]+xs_data_inc[idx]])

# 1. Total cross section interpolation
data_cs = get_data(0)
E_cs = make_axis(xs_tab_index[0][0], xs_min[0][0], xs_max[0][0], xs_spacing[0][0])
f_sigma_tot = interp1d(E_cs, data_cs, kind='cubic', bounds_error=False, fill_value=(data_cs[0], data_cs[-1]))

# 2. Scattering angle mapping R_theta = \int (1 - \cos\theta) dR
data_ang = get_data(5)
n1_ang, n2_ang = xs_tab_index[5][0], xs_tab_index[5][1]
ang_grid = data_ang.reshape(n2_ang, n1_ang)
R_rand = make_axis(n1_ang, xs_min[5][0], xs_max[5][0], xs_spacing[5][0])
E_scat = make_axis(n2_ang, xs_min[5][1], xs_max[5][1], xs_spacing[5][1])

R_theta_vals = np.zeros(n2_ang)
for i in range(n2_ang):
    theta_R = ang_grid[i, :]
    R_theta_vals[i] = integrate.simpson(1.0 - np.cos(theta_R), x=R_rand)

f_R_theta = interp1d(E_scat, R_theta_vals, kind='cubic', bounds_error=False, fill_value=(R_theta_vals[0], R_theta_vals[-1]))

# 3. I_kernel Grids
n1_I = xs_tab_index[2][0]
n2_I = xs_tab_index[2][1]
E_I = make_axis(n1_I, xs_min[2][0], xs_max[2][0], xs_spacing[2][0])
T_I = make_axis(n2_I, xs_min[2][1], xs_max[2][1], xs_spacing[2][1])

c0 = 1.3891494e6  # cm/s for 1 eV/amu

I10_new = np.zeros((n2_I, n1_I))
I11_new = np.zeros((n2_I, n1_I))
I12_new = np.zeros((n2_I, n1_I))
RR_new = np.zeros((n2_I, n1_I))

print("Calculating I_1_x and Reaction Rate integrals...", flush=True)

# Using a progress counter
total_points = n2_I * n1_I
count = 0

for iT in range(n2_I):
    for iE in range(n1_I):
        E_val = E_I[iE]
        T_val = T_I[iT]
        v_alpha = c0 * np.sqrt(E_val)
        a_beta = c0 * np.sqrt(T_val)
        delta = np.sqrt(E_val / T_val)
        
        # Integration grid bounds
        max_xi = delta + 6.0
        xi_grid = np.linspace(0, max_xi, 1000)
        
        # Relative lab energy
        Er_grid = T_val * xi_grid**2
        s1_grid = f_sigma_tot(Er_grid) * f_R_theta(Er_grid)
        
        # Define exponent terms safely
        exp_minus = np.exp(-(xi_grid - delta)**2)
        exp_plus = np.exp(-(xi_grid + delta)**2)
        
        diff_exp = exp_minus - exp_plus
        sum_exp = exp_minus + exp_plus
        
        # Integrands for I_1_x (use sigma_1 = sigma_tot * R_theta)
        integ_10 = xi_grid**2 * s1_grid * diff_exp
        integ_11 = xi_grid**3 * s1_grid * sum_exp
        integ_12 = xi_grid**4 * s1_grid * diff_exp
        
        # Integrand for Reaction Rate (uses sigma_tot only)
        # RR formula is mathematically identical to I_1_0 but with sigma_tot instead of sigma^{(1)}
        s_tot_grid = f_sigma_tot(Er_grid)
        integ_rr = xi_grid**2 * s_tot_grid * diff_exp
        
        res_10 = integrate.simpson(integ_10, x=xi_grid)
        res_11 = integrate.simpson(integ_11, x=xi_grid)
        res_12 = integrate.simpson(integ_12, x=xi_grid)
        res_rr = integrate.simpson(integ_rr, x=xi_grid)
        
        pref = 1.0 / (np.sqrt(np.pi) * v_alpha)
        I10_new[iT, iE] = pref * res_10 * (a_beta**2)
        I11_new[iT, iE] = pref * res_11 * (a_beta**3)
        I12_new[iT, iE] = pref * res_12 * (a_beta**4)
        RR_new[iT, iE]  = pref * res_rr * (a_beta**2)
        
        count += 1
        if count % 500 == 0:
            print(f"Progress: {count}/{total_points}", flush=True)

print("Integration complete. Copying back to xs_data_tab...", flush=True)

# Replace values in the xs_data_tab array
def flatten_back(arr):
    return arr.flatten()  # row-major (T then E corresponds to C-contiguous reshape which we originally achieved with reshape(n2, n1))

RR_flat  = flatten_back(RR_new)
I10_flat = flatten_back(I10_new)
I11_flat = flatten_back(I11_new)
I12_flat = flatten_back(I12_new)

idx_rr = 1; b_rr = xs_data_base[idx_rr]; inc_rr = xs_data_inc[idx_rr]
idx_10 = 2; b_10 = xs_data_base[idx_10]; inc_10 = xs_data_inc[idx_10]
idx_11 = 3; b_11 = xs_data_base[idx_11]; inc_11 = xs_data_inc[idx_11]
idx_12 = 4; b_12 = xs_data_base[idx_12]; inc_12 = xs_data_inc[idx_12]

xs_data_tab[b_rr:b_rr+inc_rr] = RR_flat
xs_data_tab[b_10:b_10+inc_10] = I10_flat
xs_data_tab[b_11:b_11+inc_11] = I11_flat
xs_data_tab[b_12:b_12+inc_12] = I12_flat

# Also we need to recalculate sigv_max since reaction_rate changed
# sigv_max is at index 6. Let's find max of RR_new
new_sigv_max = np.max(RR_new)
idx_6 = 6; b_6 = xs_data_base[idx_6]
xs_data_tab[b_6] = new_sigv_max
print(f"Updated sigv_max to {new_sigv_max:.6e}", flush=True)

print("Formatting new xs_data_tab string...", flush=True)
# To preserve formatting, we write 4 floats per line, or whatever is clean.
# Actually we can just write it as comma separated list with newlines.
formatted_strs = []
for i, val in enumerate(xs_data_tab):
    if i % 4 == 0 and i > 0:
        formatted_strs.append(f"\n    {val:.14e}")
    else:
        formatted_strs.append(f"{val:.14e}")

new_xs_tab_str = ", ".join(formatted_strs).replace("\n    , ", ",\n    ") 
# wait, better formatting:
chunks = [f"{val:.14e}" for val in xs_data_tab]
lines = []
for i in range(0, len(chunks), 4):
    lines.append("    " + ", ".join(chunks[i:i+4]))
new_xs_tab_str = ",\n".join(lines)

print("Replacing in text and writing output...", flush=True)
pattern = re.compile(r'(\bxs_data_tab\s*=\s*)(.*?)(;)', re.DOTALL)
def repl(m):
    return m.group(1) + "\n" + new_xs_tab_str + "\n  " + m.group(3)

out_text = pattern.sub(repl, orig_text)

with open(OUT_CDF, 'w') as f:
    f.write(out_text)

print(f"Successfully generated {OUT_CDF}!", flush=True)
