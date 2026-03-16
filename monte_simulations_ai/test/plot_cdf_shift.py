import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import re
import os

def load_angle_cdf(filename):
    with open(filename, 'r') as f:
        content = f.read()
    match = re.search(r'xs_data_tab\s*=\s*(.*?);', content, re.DOTALL)
    data_str = match.group(1).replace('\n', ' ').replace(',', ' ')
    tokens = data_str.split()
    data = np.array([float(x) for x in tokens])
    
    # Angle CDF is from index 10505
    angle_flat = data[10505:10505 + 51*251]
    angle_cdf = angle_flat.reshape((51, 251))
    
    energy_grid_angle = np.exp(np.linspace(np.log(0.001), np.log(100.0), 51))
    prob_grid = np.linspace(0.0, 1.0, 251)
    
    return energy_grid_angle, angle_cdf, prob_grid

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
    plt.ylabel('Scattering Angle $\chi$ (degrees)', fontsize=12)
    plt.title('Inverse Angle CDF: Probability $\\rightarrow$ Angle Mapping', fontsize=14)
    plt.grid(True, alpha=0.4)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig('inverse_cdf_comparison.png', dpi=300)
    print("Saved inverse_cdf_comparison.png")

if __name__ == '__main__':
    main()
