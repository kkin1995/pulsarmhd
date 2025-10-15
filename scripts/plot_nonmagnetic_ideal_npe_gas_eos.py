import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load data
data = pd.read_csv("/mnt/c/Users/karan/Dropbox/KARAN/2 Areas/Education/PhD/3 Research/pulsarmhd/magnetic_bps_equation_of_state.csv")
# Filter data for converged points
log_rho = data['log_rho'].values
log_P = data['log_P'].values
log_n = data['log_n'].values

# Create a combined mask to filter out negative values
mask = (log_rho > 0.0) & (log_P > 0.0) & (log_n > 0.0)

# Apply the combined mask to filter the data
log_rho = log_rho[mask]
log_P = log_P[mask]
log_n = log_n[mask]

# Calculate gamma using finite differences
def calculate_gamma_physics(log_P, log_rho):
    gamma = np.zeros_like(log_rho)

    # Use central differences for bulk
    for i in range(1, len(log_rho)-1):
        # Calculate dlog(P)/dlog(rho)
        dlogP = log_P[i+1] - log_P[i-1]
        dlogrho = log_rho[i+1] - log_rho[i-1]
        if dlogrho != 0:
            gamma[i] = dlogP/dlogrho

    # Handle endpoints using forward/backward differences
    gamma[0] = (log_P[1] - log_P[0])/(log_rho[1] - log_rho[0])
    gamma[-1] = (log_P[-1] - log_P[-2])/(log_rho[-1] - log_rho[-2])

    return gamma

# Calculate gamma
Gamma = calculate_gamma_physics(log_P, log_rho)

# Plot
plt.figure(figsize=(8, 6))
plt.scatter(log_rho, Gamma, c='green', s=10)
plt.xlabel(r'$\log_{10}(\rho)$ (g/cm$^3$)')
plt.ylabel(r'$\Gamma$')
plt.ylim(0.0, 3.0)
plt.grid(True)
plt.title('Adiabatic Index')
plt.savefig("magnetic_bps_gamma_vs_log_rho.png")

# P-rho plot
plt.figure(figsize=(8, 6))
plt.scatter(log_rho, log_P, label="B = 0 (magnetic)", color='black', s=5)
plt.xlabel(r'$\log(\rho)$ (g/cm$^3$)')
plt.ylabel(r'$\log(P)$ (dyn/cm$^2$)')
plt.grid(True)
plt.title('Magnetic BPS EOS')
plt.savefig("magnetic_bps_eos.png")
