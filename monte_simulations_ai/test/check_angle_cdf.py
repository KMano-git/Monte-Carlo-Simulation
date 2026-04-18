import numpy as np
import matplotlib.pyplot as plt
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

def main():
    cdf_path = Path(__file__).resolve().parent / '../3d3v_event-driven/run/dd_00_elastic.cdf'
    tables = load_runtime_elastic_tables(cdf_path, sigma_scale=1e-4)
    sigma_orig = tables["sigma_elastic"]
    energy_grid_sigma = tables["energy_grid_sigma"]
    angle_cdf = tables["angle_cdf"]
    energy_grid_angle = tables["energy_grid_angle"]
    prob_grid = tables["prob_grid"]
    
    def get_cx_sigma(E_eV):
        E = np.maximum(E_eV, 0.01)
        sigma = 0.6937e-18 * (1.0 - 0.155 * np.log10(E))**2 / (1.0 + 0.1112e-14 * E**3.3)
        return np.maximum(sigma, 1e-22)

    # Plot inverse CDF for a few energy points
    energies_to_plot = [0.001, 1.0, 100.0]
    
    plt.figure(figsize=(10, 6))
    for E_target in energies_to_plot:
        # Find closest energy index in angle grid
        idx = np.argmin(np.abs(energy_grid_angle - E_target))
        E_actual = energy_grid_angle[idx]
        
        # Interpolate sigma_orig for E_actual
        sig_orig = np.interp(E_actual, energy_grid_sigma, sigma_orig)
        sig_cx = get_cx_sigma(E_actual)
        
        p_cx = sig_cx / sig_orig
        
        angle_vals = np.degrees(angle_cdf[:, idx])
        
        plt.plot(prob_grid, angle_vals, label=f'{E_actual:.3f} eV (p_cx = {p_cx:.3f})')
        plt.axvline(1 - p_cx, linestyle='--', color='gray')
        
    plt.xlabel('Probability P')
    plt.ylabel('Scattering Angle $\\chi$ (degrees)')
    plt.title('Original Angle Inverse CDF')
    plt.legend()
    plt.grid()
    plt.savefig('angle_cdf_check.png')
    print('Saved angle_cdf_check.png')
    
if __name__ == '__main__':
    main()
