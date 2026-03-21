import os
import re
import numpy as np
import scipy.integrate as integrate
from scipy.interpolate import interp1d

CDF_PATH = '../3d3v_event-driven/run/dd_00_elastic.cdf'

with open(CDF_PATH, 'r') as f:
    text = f.read()

def parse_array(text, name):
    pat = re.compile(rf'\b{re.escape(name)}\s*=\s*(.*?)\s*;', re.DOTALL)
    m = pat.search(text)
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

xs_data_tab  = parse_float_array(text, 'xs_data_tab')
xs_data_base = parse_int_array(text, 'xs_data_base')
xs_data_inc  = parse_int_array(text, 'xs_data_inc')
def get_data(idx):
    return np.array(xs_data_tab[xs_data_base[idx] : xs_data_base[idx]+xs_data_inc[idx]])

min_f = parse_float_array(text, 'xs_min'); xs_min = [min_f[i*3:(i+1)*3] for i in range(20)]
spc_r = parse_array(text, 'xs_spacing'); xs_spacing = [spc_r[i*4:(i+1)*4] for i in range(20)]
tab_f = parse_int_array(text, 'xs_tab_index'); xs_tab_index = [tab_f[i*3:(i+1)*3] for i in range(20)]
max_f = parse_float_array(text, 'xs_max'); xs_max = [max_f[i*3:(i+1)*3] for i in range(20)]

data_cs = get_data(0)
E_cs = make_axis(xs_tab_index[0][0], xs_min[0][0], xs_max[0][0], xs_spacing[0][0])
f_sigma_tot = interp1d(E_cs, data_cs, kind='cubic', bounds_error=False, fill_value=(data_cs[0], data_cs[-1]))

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

data_I10 = get_data(2)
I10_ex = data_I10.reshape(xs_tab_index[2][1], xs_tab_index[2][0])
E_I = make_axis(xs_tab_index[2][0], xs_min[2][0], xs_max[2][0], xs_spacing[2][0])
T_I = make_axis(xs_tab_index[2][1], xs_min[2][1], xs_max[2][1], xs_spacing[2][1])

data_rr = get_data(1).reshape(xs_tab_index[1][1], xs_tab_index[1][0])

c0 = 1.3891494e6

def test_er_scaling(scale):
    val_calc = []
    val_exp = []
    
    for iT in [10, 25, 40]:
        for iE in [10, 25, 40]:
            E_val = E_I[iE]
            T_val = T_I[iT]
            v_alpha = c0 * np.sqrt(E_val)
            a_beta = c0 * np.sqrt(T_val)
            delta = np.sqrt(E_val / T_val)

            max_xi = delta + 6.0
            xi_grid = np.linspace(0, max_xi, 1000)
            Er_grid = scale * T_val * xi_grid**2
            s1_grid = f_sigma_tot(Er_grid) * f_R_theta(Er_grid)
            integ_grid = xi_grid**2 * s1_grid * (np.exp(-(xi_grid - delta)**2) - np.exp(-(xi_grid + delta)**2))
            res = integrate.simpson(integ_grid, x=xi_grid)
            calc = (a_beta**2 / (np.sqrt(np.pi) * v_alpha)) * res
            exp = I10_ex[iT, iE]
            
            # reaction rate
            s_tot_grid = f_sigma_tot(Er_grid)
            rr_integ = xi_grid**2 * s_tot_grid * (np.exp(-(xi_grid - delta)**2) - np.exp(-(xi_grid + delta)**2))
            rr_res = integrate.simpson(rr_integ, x=xi_grid)
            calc_rr = (a_beta**2 / (np.sqrt(np.pi) * v_alpha)) * rr_res
            exp_rr = data_rr[iT, iE]

            print(f"E={E_val:8.2e}, T={T_val:8.2e} | I10 Ratio: {calc/exp:6.3f} | RR Ratio: {calc_rr/exp_rr:6.3f}")  

test_er_scaling(1.0)

# Also check sigma_1 scaling?
# Or maybe the integration formula has a different prefactor?
