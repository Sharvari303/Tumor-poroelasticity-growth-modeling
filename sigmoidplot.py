import numpy as np
import matplotlib.pyplot as plt

# Hill function parameters
hill_coefficient = 2.0
K_half = 1e-8  # Half-max concentration (M)
max_prob = 2e-4  # Saturation level

# Growth factor concentration range (log scale)
gf_conc = np.logspace(-9.5, -7, 1000)

# Hill function
def hill_function(conc, K_half, n, max_prob):
    return max_prob * (conc**n) / (K_half**n + conc**n)

# Compute probability values
prob = hill_function(gf_conc, K_half, hill_coefficient, max_prob)

# Plot
plt.figure(figsize=(10, 8))
plt.plot(gf_conc, prob, label=f'Hill Function (Hill Coeff = {hill_coefficient})', color='blue')

# Add horizontal saturation line
plt.axhline(y=max_prob, color='green', linestyle='--', label='Saturation (2e-4)')

# Add vertical line at half-max concentration
plt.axvline(x=K_half, color='red', linestyle='--', label='Half-max (1e-8)')

# Plot formatting
plt.xscale('log')
plt.xlabel('Growth Factor Concentration (M)')
plt.ylabel('Cellular Growth Rate')
plt.title(f'ABM cell growth rate - sigmoidal function of GF concentration')
plt.legend()
plt.grid(True)
plt.grid(which='major', linestyle='-', linewidth=0.5, alpha=0.7)
plt.minorticks_on()
plt.grid(which='minor', axis = 'x',linestyle='-', linewidth=0.5, alpha=0.4)

plt.savefig('sigmoidplot.png')



