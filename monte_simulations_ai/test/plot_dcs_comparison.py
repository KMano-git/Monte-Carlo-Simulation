import csv
import matplotlib.pyplot as plt
import numpy as np

# Load DCS data from test_dd_00 output (calculated utilizing Fortran logic against pure EL CDF)
energies = []
data = {}

with open('test_dd_00/dcs_output.csv', 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        energy = float(row['Energy_eV'])
        angle = float(row['Angle_deg'])
        dcs = float(row['DCS_m2_sr'])
        
        if energy not in data:
            data[energy] = {'angles': [], 'dcs': []}
            if energy not in energies:
                energies.append(energy)
            
        data[energy]['angles'].append(angle)
        data[energy]['dcs'].append(dcs)


# Functions to calculate theoretical cross-sections
def get_cx_sigma(E_eV):
    """Calculate the total charge exchange cross-section [m^2] according to Janev's formula"""
    E = max(float(E_eV), 0.01)
    sigma = 0.6937e-18 * (1.0 - 0.155 * np.log10(E))**2 / (1.0 + 0.1112e-14 * E**3.3)
    return max(sigma, 1e-22)

# Create theoretical DCS curve for CX
def janev_dcs(angle_deg, energy_eV):
    """
    Theoretical DCS for Charge Exchange.
    Usually we assign backwards scattering for charge exchange in the CM frame.
    Here we plot a strong backscatter peaked distribution that integrates to total_sigma.
    Model: DCS ~ exp(-\theta_diff) where \theta_diff = 180 - \theta
    """
    total_sigma = get_cx_sigma(energy_eV)
    
    # Simple exponential model peaked at 180 degrees
    # e^(alpha * (cos(theta) - 1)) / normalization, but peaked at backward so cos(theta) + 1
    # Let's use alpha = 10 for a sharp peak
    angle_rad = np.radians(angle_deg)
    alpha = 10.0
    
    # Normalization for e^(alpha * cos(theta)) from -1 to 1 is 
    # 2*pi * integral_(-1)^1 e^(alpha*x) dx = 2*pi * (e^alpha - e^-alpha) / alpha
    norm = 2 * np.pi * (np.exp(alpha) - np.exp(-alpha)) / alpha
    
    # We want it peaked at 180 (cos(180) = -1), so we use -alpha * cos(theta)
    dcs = (total_sigma / norm) * np.exp(-alpha * np.cos(angle_rad))
    return dcs

plt.figure(figsize=(12, 8))

colors = ['b', 'g', 'r', 'c', 'm']
for i, energy in enumerate(energies):
    color = colors[i % len(colors)]
    
    # Plotting modified EL DCS
    plt.plot(data[energy]['angles'], data[energy]['dcs'], 
             label=f'Modified Elastic ({energy} eV)', color=color)
    
    # Plotting Janev CX DCS
    angles = np.linspace(0.1, 179.9, 100)
    cx_dcs_vals = [janev_dcs(a, energy) for a in angles]
    
    plt.plot(angles, cx_dcs_vals, 
             label=f'Janev CX (Backscatter Model) ({energy} eV)', color=color, linestyle='--')

plt.yscale('log')
plt.xlabel('Scattering Angle (degrees)')
plt.ylabel('Differential Cross Section (m^2/sr)')
plt.title('Comparison of Modified Elastic DCS vs Janev CX DCS (Strong Backscatter peaked)')
plt.grid(True, which='both', ls='--', alpha=0.5)
plt.xlim(0, 180)
plt.legend(bbox_to_anchor=(1.01, 1), loc='upper left')
plt.tight_layout()

plt.savefig('comparative_dcs_plot.png', dpi=300)
print('Plot saved to comparative_dcs_plot.png')
