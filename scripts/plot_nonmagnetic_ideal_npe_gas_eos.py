import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load data
data = pd.read_csv("/mnt/c/Users/karan/Dropbox/KARAN/2 Areas/Education/PhD/3 Research/pulsarmhd/ideal_non_magnetic_npe_gas_eos.csv")
log_rho = data['log_rho']
log_P = data['log_P']
log_n = data['log_n']
converged = data['converged']

Gamma = np.gradient(log_P, log_n)

# Plot log_P vs log_rho
plt.figure(figsize=(8, 6))
plt.scatter(log_rho, log_P, label="B = 0 (nonmagnetic)", color='black', s=5)
plt.xlabel(r'$\log(\rho)$ (g/cm$^3$)')
plt.ylabel(r'$\log(P)$ (dyn/cm$^2$)')
plt.grid(True)
plt.legend()
plt.savefig("ideal_non_magnetic_npe_gas_eos.png")

# # Plot Gamma vs log_rho
# plt.figure(figsize=(8, 6))
# plt.scatter(log_rho, Gamma, label=r'$\Gamma$', color='blue', s=5)
# plt.xlabel(r'$\log(\rho)$ (g/cm$^3$)')
# plt.ylabel(r'$\Gamma$')
# plt.grid(True)
# plt.legend()
# plt.savefig("gamma_vs_log_rho.png")

plt.figure(figsize=(8, 6))
plt.scatter(log_rho[converged == 1], Gamma[converged == 1], c='green', label='Converged', s=10)
plt.scatter(log_rho[converged == 0], Gamma[converged == 0], c='red', label='Failed', s=10)
plt.xlabel(r'$\log_{10}(\rho)$ (g/cm$^3$)')
plt.ylabel(r'$\Gamma$')
plt.legend()
plt.grid(True)
plt.title('Adiabatic Index')
plt.savefig("gamma_vs_log_rho.png")
