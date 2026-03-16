import numpy as np
import matplotlib.pyplot as plt
import re

def main():
    with open('../3d3v_event-driven/run/dd_00_elastic.cdf', 'r') as f:
        content = f.read()
        
    match = re.search(r'xs_data_tab\s*=\s*(.*?);', content, re.DOTALL)
    data_str = match.group(1)
    tokens = re.findall(r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?', data_str)
    data = np.array([float(x) for x in tokens])
    
    # 1. 101 points sigma_orig (cm^2)
    sigma_orig = data[:101] * 1e-4 # m^2
    energy_grid_sigma = np.logspace(np.log10(0.001), np.log10(100), 101)
    
    # 6. 251 x 51 angle CDF
    angle_cdf_flat = data[10505:10505 + 251 * 51]
    angle_cdf = angle_cdf_flat.reshape((51, 251)).T # Fortran is column major: (251, 51)
    
    energy_grid_angle = np.logspace(np.log10(0.001), np.log10(100), 51)
    prob_grid = np.linspace(0, 1, 251)
    
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
    plt.ylabel('Scattering Angle \chi (degrees)')
    plt.title('Original Angle Inverse CDF')
    plt.legend()
    plt.grid()
    plt.savefig('angle_cdf_check.png')
    print('Saved angle_cdf_check.png')
    
if __name__ == '__main__':
    main()
