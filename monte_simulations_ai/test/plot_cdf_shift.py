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

def load_angle_cdf(filename):
    cdf_path = Path(filename)
    if not cdf_path.exists():
        cdf_path = Path(__file__).resolve().parent / filename
    tables = load_runtime_elastic_tables(cdf_path)
    return tables["energy_grid_angle"], tables["angle_cdf"].T, tables["prob_grid"]

def main():
    orig_file = '../3d3v_event-driven/run/dd_00_elastic.cdf'
    new_file = '../3d3v_event-driven/run/dd_00_elastic_pure_el_angle.cdf'
    
    E_e_ang, A_orig, P_grid = load_angle_cdf(orig_file)
    _, A_new, _ = load_angle_cdf(new_file)
    
    energies_to_plot = [0.1, 1.0, 10.0, 100.0]
    
    # Plot 1: Inverse CDF (Angle vs Probability)
    plt.figure(figsize=(12, 6))
    
    colors = ['b', 'g', 'r', 'c']
    for i, e_target in enumerate(energies_to_plot):
        idx = np.argmin(np.abs(E_e_ang - e_target))
        actual_e = E_e_ang[idx]
        
        c = colors[i]
        
        # Original (dashed)
        angles_orig_deg = np.degrees(A_orig[idx, :])
        plt.plot(P_grid, angles_orig_deg, c=c, linestyle='--', alpha=0.7, 
                 label=f'Original ({actual_e:.1f} eV)')
        
        # New Pure EL (solid)
        angles_new_deg = np.degrees(A_new[idx, :])
        plt.plot(P_grid, angles_new_deg, c=c, linestyle='-', linewidth=2,
                 label=f'Pure EL ({actual_e:.1f} eV)')
                 
    plt.xlabel('Cumulative Probability $P$', fontsize=12)
    plt.ylabel('Scattering Angle $\\chi$ (degrees)', fontsize=12)
    plt.title('Inverse Angle CDF: Probability $\\rightarrow$ Angle Mapping', fontsize=14)
    plt.grid(True, alpha=0.4)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig('inverse_cdf_comparison.png', dpi=300)
    print("Saved inverse_cdf_comparison.png")

if __name__ == '__main__':
    main()
