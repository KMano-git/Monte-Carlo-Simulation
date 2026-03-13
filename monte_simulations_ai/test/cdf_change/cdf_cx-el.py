import numpy as np
import os
import re

def sigma_cx(E_eV):
    """
    Calculate Charge Exchange (CX) cross-section using Janev's formula.
    Returns the cross-section in cm^2.
    """
    E_safe = max(float(E_eV), 0.01)
    
    # Janev formula [m^2]
    sigma_m2 = 0.6937e-18 * (1.0 - 0.155 * np.log10(E_safe))**2 / (1.0 + 0.1112e-14 * E_safe**3.3)
    
    # Apply minimum bound
    sigma_m2 = max(sigma_m2, 1.0e-22)
    
    # Convert m^2 to cm^2
    sigma_cm2 = sigma_m2 * 1e4
    return sigma_cm2

def main():
    # File paths relative to this script's location (test/cdf_change/)
    input_cdf = 'dd_00_elastic.cdf'
    output_cdf = 'dd_00_elastic_pure_el.cdf'
    
    if not os.path.exists(input_cdf):
        print(f"Error: Inut file {input_cdf} not found.")
        return
        
    print(f"Reading original CDF file: {input_cdf}")
    with open(input_cdf, 'r') as f:
        content = f.read()
        
    # Extract the numerical data enclosed between xs_data_tab = ... ;
    match = re.search(r'xs_data_tab\s*=\s*(.*?);', content, re.DOTALL)
    if not match:
        print("Error: 'xs_data_tab = ... ;' section not found in the CDF file.")
        return
        
    data_str = match.group(1)
    
    # Extract all floating point numbers using regular expressions
    # This handles integers, floats, and scientific notation
    tokens = re.findall(r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?', data_str)
    
    if len(tokens) < 23320:
        print(f"Error: Found only {len(tokens)} values, expected 23320.")
        return
        
    # Energy grid setup (same as Fortran code)
    N_ENERGY_SIGMA = 101
    log_E_min = np.log(0.001)
    log_E_max = np.log(100.0)
    log_E_grid_dlt = (log_E_max - log_E_min) / (N_ENERGY_SIGMA - 1.0)
    
    modified_data = []
    
    print("\nModifying cross sections (Original EL -> Pure EL)...")
    print(f"{'Energy_eV':>10} | {'Orig_EL_cm2':>15} | {'CX_cm2':>15} | {'New_EL_cm2':>15}")
    print("-" * 65)
    
    for i in range(23320):
        val = float(tokens[i])
        
        if i < N_ENERGY_SIGMA:
            # Calculate energy [eV] at this index
            E_grid = np.exp(log_E_min + i * log_E_grid_dlt)
            
            # Calculate CX cross-section [cm^2]
            cx_cm2 = sigma_cx(E_grid)
            
            # Original subtraction
            raw_new_val = val - cx_cm2
            
            # Check for negative values
            if raw_new_val < 0.0:
                print(f"WARNING: Negative cross-section detected at Energy = {E_grid:.5f} eV.")
                print(f"         Orig EL = {val:.5e}, CX = {cx_cm2:.5e}, Diff = {raw_new_val:.5e}")
                print(f"         Setting New EL to 0.0")
                new_val = 0.0
            else:
                new_val = raw_new_val
            
            if i % 10 == 0:
                print(f"{E_grid:10.5f} | {val:15.5e} | {cx_cm2:15.5e} | {new_val:15.5e}")
                
            modified_data.append(new_val)
        else:
            modified_data.append(val)
            
    # Format the updated numerical array into chunks of 5 per line
    lines = []
    for i in range(0, 23320, 5):
        chunk = modified_data[i:i+5]
        # Format values in scientific notation with high precision
        line = "    " + ", ".join(f"{x:.14e}" for x in chunk)
        lines.append(line)
        
    # Reassemble the string
    new_data_str = ",\n".join(lines)
    
    # Inject it back into the content (keeping the rest of the NetCDF structure intact)
    new_content = content[:match.start()] + "xs_data_tab =\n" + new_data_str + "\n;" + content[match.end():]
    
    # Write the modified data back to a new file
    with open(output_cdf, 'w') as f:
        f.write(new_content)
        
    print(f"\nSuccessfully wrote pure EL modified CDF to: {output_cdf}")

if __name__ == '__main__':
    main()
