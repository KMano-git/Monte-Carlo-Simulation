#!/usr/bin/env python3
"""Plot solution of 1D Poisson equation"""

import numpy as np
import matplotlib.pyplot as plt

# Load data
data = np.loadtxt('poisson_result.dat', delimiter=',')
x = data[:, 0]
V = data[:, 1]

# Plot settings
plt.figure(figsize=(10, 6))
plt.plot(x, V, 'b-', linewidth=2, marker='o', markersize=4, label='V(x)')

plt.xlabel('x', fontsize=14)
plt.ylabel('V (Voltage)', fontsize=14)
plt.title('1D Poisson Equation Solution', fontsize=16)
plt.grid(True, alpha=0.3)
plt.legend(fontsize=12)

# Save and show
plt.tight_layout()
plt.savefig('poisson_plot.png', dpi=150)
plt.show()
print('Plot saved to poisson_plot.png')
