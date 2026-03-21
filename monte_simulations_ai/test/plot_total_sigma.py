import numpy as np
import matplotlib.pyplot as plt
import re
import os

def get_cx_sigma(E_eV):
    """
    Calculate Charge Exchange (CX) cross-section using Janev's formula.
    Returns the cross-section in m^2.
    """
    E = np.maximum(E_eV, 0.01)
    sigma = 0.6937e-18 * (1.0 - 0.155 * np.log10(E))**2 / (1.0 + 0.1112e-14 * E**3.3)
    return np.maximum(sigma, 1e-22)

def load_cdf_sigma(filename):
    """Load the first 101 values from the CDF xs_data_tab (which are the cross sections)"""
    with open(filename, 'r') as f:
        content = f.read()
    
    match = re.search(r'xs_data_tab\s*=\s*(.*?);', content, re.DOTALL)
    if not match:
        raise ValueError(f"xs_data_tab not found in {filename}")
        
    data_str = match.group(1)
    tokens = re.findall(r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?', data_str)
    
    # Return the first 101 items (convert cm^2 to m^2)
    # The pure CDF and original CDF have cm^2 data
    return np.array([float(x) for x in tokens[:101]]) * 1e-4

def load_cdf_energy_grid(filename):
    """Load the energy grid from CDF xs_min/xs_max for cross_section (index 0)"""
    with open(filename, 'r') as f:
        content = f.read()
    
    # xs_min は (20, 3) の配列: 最初の3要素が cross_section の min
    min_match = re.search(r'xs_min\s*=\s*(.*?);', content, re.DOTALL)
    min_vals = re.findall(r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?', min_match.group(1))
    E_min = float(min_vals[0])  # cross_section の 1st axis min
    
    max_match = re.search(r'xs_max\s*=\s*(.*?);', content, re.DOTALL)
    max_vals = re.findall(r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?', max_match.group(1))
    E_max = float(max_vals[0])  # cross_section の 1st axis max
    
    return np.logspace(np.log10(E_min), np.log10(E_max), 101)

def main():
    # File paths
    orig_cdf_path = '../3d3v_event-driven/run/dd_00_elastic.cdf'
    pure_cdf_path = '../3d3v_event-driven/run/dd_00_elastic_pure_el.cdf'
    
    # Energy grid をCDFメタデータから読む (eV/amu 単位)
    energy_grid = load_cdf_energy_grid(orig_cdf_path)
    
    # 1. Load Original EL Cross-section (m^2)
    sigma_orig_el = load_cdf_sigma(orig_cdf_path)
    
    # 2. Load Pure EL Cross-section (m^2)
    sigma_pure_el = load_cdf_sigma(pure_cdf_path)
    
    # 3. Calculate Janev CX Cross-section (m^2)
    # CX断面積の引数はeV単位: D+D⁺系では eV/amu ≈ eV (μ ≈ 1 amu)
    sigma_cx = get_cx_sigma(energy_grid)
    
    # Plotting
    plt.figure(figsize=(10, 7))
    
    plt.plot(energy_grid, sigma_orig_el, 'k-', linewidth=3, label='Original CDF (Experiment / EL + CX)')
    plt.plot(energy_grid, sigma_pure_el, 'b-', linewidth=2, label='Modified Pure Elastic (EL)')
    plt.plot(energy_grid, sigma_cx, 'r--', linewidth=2, label='Janev Formula (CX)')
    
    # We can also plot the sum of Pure EL + CX just to verify it matches original
    plt.plot(energy_grid, sigma_pure_el + sigma_cx, 'y:', linewidth=2, label='Sum: Pure EL + CX')
    
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(energy_grid[0], energy_grid[-1])
    plt.ylim(1e-19, 1e-16)
    
    plt.xlabel('Energy (eV/amu)', fontsize=14)
    plt.ylabel('Total Cross-Section (m$^2$)', fontsize=14)
    plt.title('Total Scattering Cross-Sections: Elastic vs Charge Exchange', fontsize=15)
    plt.grid(True, which="both", ls="-", alpha=0.2)
    plt.legend(fontsize=12)
    
    plt.tight_layout()
    output_filename = 'total_cross_section_comparison.png'
    plt.savefig(output_filename, dpi=300)
    print(f'Plot saved successfully to {output_filename}')

if __name__ == '__main__':
    main()
