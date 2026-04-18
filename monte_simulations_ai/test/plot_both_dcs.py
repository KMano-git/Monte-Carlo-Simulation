import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
from pathlib import Path
import sys


def _resolve_el_data_dir():
    for parent in Path(__file__).resolve().parents:
        candidate = parent / "el_data" / "cdf_compat.py"
        if candidate.exists():
            return candidate.parent
    raise FileNotFoundError("Could not locate el_data/cdf_compat.py")


EL_DATA_DIR = _resolve_el_data_dir()
if str(EL_DATA_DIR) not in sys.path:
    sys.path.insert(0, str(EL_DATA_DIR))

from cdf_compat import load_runtime_elastic_tables

def load_cdf(filename):
    cdf_path = Path(filename)
    if not cdf_path.exists():
        cdf_path = Path(__file__).resolve().parent / filename
    tables = load_runtime_elastic_tables(cdf_path, sigma_scale=1e-4)
    return (
        tables["energy_grid_sigma"],
        tables["sigma_elastic"],
        tables["energy_grid_angle"],
        tables["angle_cdf"].T,
        tables["prob_grid"],
    )

def get_dcs(energy, energy_grid_sigma, sigma, energy_grid_angle, angle_cdf, prob_grid):
    sig_val = np.interp(energy, energy_grid_sigma, sigma)
    
    # Interpolate angle CDF for this energy
    idx = np.searchsorted(energy_grid_angle, energy)
    if idx == 0:
        angle_vals = angle_cdf[0, :]
    elif idx == len(energy_grid_angle):
        angle_vals = angle_cdf[-1, :]
    else:
        # Linear interp on log energy
        t = (np.log(energy) - np.log(energy_grid_angle[idx-1])) / \
            (np.log(energy_grid_angle[idx]) - np.log(energy_grid_angle[idx-1]))
        angle_vals = (1-t) * angle_cdf[idx-1, :] + t * angle_cdf[idx, :]
        
    angles_mid = []
    dcs_vals = []
    for i in range(1, len(prob_grid)):
        p1 = prob_grid[i-1]
        p2 = prob_grid[i]
        
        p1 = max(1e-6, p1)
        p2 = min(1.0 - 1e-6, p2)
        
        chi1 = np.interp(p1, prob_grid, angle_vals)
        chi2 = np.interp(p2, prob_grid, angle_vals)
        
        dchi = chi2 - chi1
        if dchi > 1e-12:
            chi_mid = (chi1 + chi2) / 2.0
            if np.sin(chi_mid) > 1e-12:
                dP_dchi = (p2 - p1) / dchi
                dcs = (sig_val / (2.0 * np.pi * np.sin(chi_mid))) * dP_dchi
                angles_mid.append(np.degrees(chi_mid))
                dcs_vals.append(dcs)
                
    return np.array(angles_mid), np.array(dcs_vals)

def main():
    orig_file = '../3d3v_event-driven/run/dd_00_elastic.cdf'
    new_file = '../3d3v_event-driven/run/dd_00_elastic_pure_el_angle.cdf'
    
    E_e_sig, S_orig, E_e_ang, A_orig, P_grid = load_cdf(orig_file)
    E_e_sig2, S_new, E_e_ang2, A_new, P_grid2 = load_cdf(new_file)
    
    energies = [0.1, 1.0, 10.0, 100.0]
    
    plt.figure(figsize=(10, 6))
    colors = ['b', 'g', 'r', 'c']
    for i, e in enumerate(energies):
        c = colors[i]
        
        # Original
        ang_orig, dcs_orig = get_dcs(e, E_e_sig, S_orig, E_e_ang, A_orig, P_grid)
        plt.plot(ang_orig, dcs_orig, c=c, linestyle='-', label=f'Orig {e} eV')
        
        # New
        ang_new, dcs_new = get_dcs(e, E_e_sig2, S_new, E_e_ang2, A_new, P_grid2)
        plt.plot(ang_new, dcs_new, c=c, linestyle='--', linewidth=3, alpha=0.7, label=f'New Pure {e} eV')
        
    plt.yscale('log')
    plt.xlim(0, 180)
    plt.xlabel('Angle (deg)')
    plt.ylabel('DCS (m^2/sr)')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, which='both', alpha=0.3)
    plt.tight_layout()
    plt.savefig('true_dcs_comparison.png', dpi=300)
    print("Saved true_dcs_comparison.png")

if __name__ == "__main__":
    main()
