import numpy as np
import re
import os

def sigma_cx(E_eV):
    """
    Calculate Charge Exchange (CX) cross-section using Janev's formula.
    Returns the cross-section in cm^2.
    """
    E_safe = np.maximum(np.array(E_eV, dtype=float), 0.01)
    
    # Janev formula [m^2]
    sigma_m2 = 0.6937e-18 * (1.0 - 0.155 * np.log10(E_safe))**2 / (1.0 + 0.1112e-14 * E_safe**3.3)
    sigma_m2 = np.maximum(sigma_m2, 1.0e-22)
    
    # Convert m^2 to cm^2
    return sigma_m2 * 1e4

def main():
    input_cdf = '../../3d3v_event-driven/run/dd_00_elastic.cdf'
    output_cdf = '../../3d3v_event-driven/run/dd_00_elastic_pure_el_angle.cdf'
    
    if not os.path.exists(input_cdf):
        print(f"Error: Input file {input_cdf} not found.")
        return
        
    print(f"Reading original CDF file: {input_cdf}")
    with open(input_cdf, 'r') as f:
        content = f.read()
        
    match = re.search(r'xs_data_tab\s*=\s*(.*?);', content, re.DOTALL)
    if not match:
        print("Error: 'xs_data_tab = ... ;' section not found.")
        return
        
    data_str = match.group(1)
    tokens = re.findall(r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?', data_str)
    
    if len(tokens) < 23320:
        print(f"Error: Found only {len(tokens)} values.")
        return
        
    data = np.array([float(x) for x in tokens])
    
    # 1. Update total cross-sections (Elements 0 to 100)
    N_ENERGY_SIGMA = 101
    log_E_sigma_min = np.log(0.001)
    log_E_sigma_max = np.log(100.0)
    E_grid_sigma = np.exp(np.linspace(log_E_sigma_min, log_E_sigma_max, N_ENERGY_SIGMA))
    
    sigma_orig = data[0:N_ENERGY_SIGMA].copy()
    cx_cm2_sigma = sigma_cx(E_grid_sigma)
    sigma_pure_el = np.maximum(0.0, sigma_orig - cx_cm2_sigma)
    
    # 2. Update Angle CDF (Elements 10505 to 23305)
    N_ENERGY_ANGLE = 51
    N_PROB = 251
    
    log_E_angle_min = np.log(0.001)
    log_E_angle_max = np.log(100.0)
    E_grid_angle = np.exp(np.linspace(log_E_angle_min, log_E_angle_max, N_ENERGY_ANGLE))
    
    prob_grid = np.linspace(0.0, 1.0, N_PROB)
    
    angle_cdf_orig_flat = data[10505:10505 + N_ENERGY_ANGLE * N_PROB]
    # Reshape to (N_ENERGY_ANGLE, N_PROB) -> Each row is one energy
    angle_cdf_orig = angle_cdf_orig_flat.reshape((N_ENERGY_ANGLE, N_PROB))
    angle_cdf_new  = np.zeros_like(angle_cdf_orig)
    
    print("\nModifying Angle CDFs...")
    for j in range(N_ENERGY_ANGLE):
        E_j = E_grid_angle[j]
        # Interpolate original sigma at E_j
        sig_orig_j = np.interp(E_j, E_grid_sigma, sigma_orig)
        sig_cx_j = sigma_cx(E_j)
        
        p_cx = sig_cx_j / sig_orig_j if sig_orig_j > 1e-30 else 0.0
        p_cx = np.clip(p_cx, 0.0, 0.9999) # avoid division by zero
        
        if j % 10 == 0:
            print(f"Energy: {E_j:8.3f} eV, Orig: {sig_orig_j:.3e}, CX: {sig_cx_j:.3e}, p_cx: {p_cx*100:.1f} %")
            
        # Since the angle CDF in Fortran decreases from 180 deg at P=0 to 0 deg at P=1
        # The CX component (180 deg) is concentrated at P in [0, p_cx].
        # We need to preserve the elastic component which is at P in [p_cx, 1.0].
        # So for P_new in [0, 1], P_orig = p_cx + P_new * (1.0 - p_cx)
        P_orig_target = p_cx + prob_grid * (1.0 - p_cx)
        
        # Interpolate the angles corresponding to these probabilities
        angles_orig_j = angle_cdf_orig[j, :]
        angle_cdf_new[j, :] = np.interp(P_orig_target, prob_grid, angles_orig_j)
    
    # 3. Reconstruct output data array
    modified_data = data.copy()
    modified_data[0:N_ENERGY_SIGMA] = sigma_pure_el
    modified_data[10505:10505 + N_ENERGY_ANGLE * N_PROB] = angle_cdf_new.flatten()
    
    # 4. Format and save
    lines = []
    for i in range(0, 23320, 5):
        chunk = modified_data[i:i+5]
        line = "    " + ", ".join(f"{x:.14e}" for x in chunk)
        lines.append(line)
        
    new_data_str = ",\n".join(lines)
    new_content = content[:match.start()] + "xs_data_tab =\n" + new_data_str + "\n;" + content[match.end():]
    
    with open(output_cdf, 'w') as f:
        f.write(new_content)
        
    print(f"\nSuccessfully wrote pure EL modified CDF (with Angle correction) to: {output_cdf}")

if __name__ == '__main__':
    main()
