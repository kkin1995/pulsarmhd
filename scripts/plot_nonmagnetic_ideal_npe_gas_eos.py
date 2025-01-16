import pandas as pd
import matplotlib.pyplot as plt

# Load data
data = pd.read_csv("/mnt/c/Users/karan/Dropbox/KARAN/2 Areas/Education/PhD/3 Research/pulsarmhd/ideal_non_magnetic_npe_gas_eos.csv")
log_rho = data['log_rho']
log_P = data['log_P']

# Plot
plt.figure(figsize=(8, 6))
plt.scatter(log_rho, log_P, label="B = 0 (nonmagnetic)", color='black', s = 5)

# Labels and title
plt.xlabel(r'$\log(\rho)$ (g/cm$^3$)')
plt.ylabel(r'$\log(P)$ (dyn/cm$^2$)')
plt.grid(True)
plt.legend()
plt.savefig("ideal_non_magnetic_npe_gas_eos.png")
