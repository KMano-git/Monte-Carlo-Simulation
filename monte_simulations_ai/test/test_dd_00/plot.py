import csv
import matplotlib.pyplot as plt

energies = []
data = {}

with open('dcs_output.csv', 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        energy = float(row['Energy_eV'])
        angle = float(row['Angle_deg'])
        dcs = float(row['DCS_m2_sr'])
        
        if energy not in data:
            data[energy] = {'angles': [], 'dcs': []}
            energies.append(energy)
            
        data[energy]['angles'].append(angle)
        data[energy]['dcs'].append(dcs)

plt.figure(figsize=(10, 6))

for energy in energies:
    plt.plot(data[energy]['angles'], data[energy]['dcs'], label=f'{energy} eV')

plt.yscale('log')
plt.xlabel('Scattering Angle (degrees)')
plt.ylabel('Differential Cross Section (m^2/sr)')
plt.title('Elastic Scattering Differential Cross Section')
plt.grid(True, which='both', ls='--', alpha=0.5)
plt.legend()
plt.tight_layout()
plt.savefig('dcs_plot.png', dpi=300)
print('Plot saved to dcs_plot.png')
