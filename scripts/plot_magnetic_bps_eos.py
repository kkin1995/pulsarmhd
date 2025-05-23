import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema

# Define a function to calculate Gamma using finite differences
def calculate_gamma_physics(log_P, log_rho):
    gamma = np.zeros_like(log_rho)
    
    # Use central differences for interior points
    for i in range(1, len(log_rho)-1):
        dlogP = log_P[i+1] - log_P[i-1]
        dlogrho = log_rho[i+1] - log_rho[i-1]
        if dlogrho != 0:
            gamma[i] = dlogP / dlogrho
    
    # Use forward and backward differences at the endpoints
    gamma[0] = (log_P[1] - log_P[0]) / (log_rho[1] - log_rho[0])
    gamma[-1] = (log_P[-1] - log_P[-2]) / (log_rho[-1] - log_rho[-2])
    
    return gamma

# List of file suffixes corresponding to different magnetic field values
file_suffixes = ['001', '01', '1', '10', '100']
# Create a dictionary mapping file suffix to magnetic field value (for labeling)
magnetic_fields = {
    '001': '0.01',
    '01': '0.1',
    '1': '1.0',
    '10': '10.0',
    '100': '100.0'
}

# Prepare figures for EOS and Gamma plots
plt.figure(figsize=(8, 6))
plt.title("Magnetic BPS EOS (log(P) vs. log(ρ))")
plt.xlabel(r'$\log_{10}(\rho)$ (g/cm$^3$)')
plt.ylabel(r'$\log_{10}(P)$ (dyn/cm$^2$)')
plt.grid(True)

plt.figure(figsize=(8, 6))
plt.title("Adiabatic Index (Γ) vs. log(ρ)")
plt.xlabel(r'$\log_{10}(\rho)$ (g/cm$^3$)')
plt.ylabel(r'$\Gamma$')
plt.ylim(0.0, 3.0)
plt.grid(True)

# We will use separate figures but later display/save them separately.
# For clarity, we will create two figures and add each dataset to both.
fig_eos = plt.figure("EOS", figsize=(8, 6))
ax_eos = fig_eos.add_subplot(1,1,1)
ax_eos.set_xlabel(r'$\log_{10}(\rho)$ (g/cm$^3$)')
ax_eos.set_ylabel(r'$\log_{10}(P)$ (dyn/cm$^2$)')
ax_eos.set_title("Magnetic BPS EOS")
ax_eos.grid(True)

fig_gamma = plt.figure("Gamma", figsize=(8, 6))
ax_gamma = fig_gamma.add_subplot(1,1,1)
ax_gamma.set_xlabel(r'$\log_{10}(\rho)$ (g/cm$^3$)')
ax_gamma.set_ylabel(r'$\Gamma$')
ax_gamma.set_title("Adiabatic Index (Γ) vs. log(ρ)")
ax_gamma.set_ylim(0.0, 5.0)
ax_gamma.grid(True)

# Define a list of colors for plotting
colors = ['black', 'red', 'blue', 'green', 'orange']

# Loop over each file suffix
for idx, suffix in enumerate(file_suffixes):
    filename = f"magnetic_bps_eos_B_{suffix}.csv"
    
    # Load data from the CSV file
    data = pd.read_csv(filename)
    
    # Extract columns
    log_rho = data['log_rho'].values
    log_P = data['log_P'].values
    log_n = data['log_n'].values

    # Create mask to filter out any negative values
    mask = (log_rho > 0.0) & (log_P > 0.0) & (log_n > 0.0)
    log_rho = log_rho[mask]
    log_P = log_P[mask]
    log_n = log_n[mask]
    
    # Calculate the adiabatic index Gamma
    Gamma = calculate_gamma_physics(log_P, log_rho)
    
    # Plot log(P) vs. log(ρ) on the EOS figure
    ax_eos.scatter(log_rho, log_P, s=5, color=colors[idx],
                   label=f"B = {magnetic_fields[suffix]}")
    
    # Plot only the upper envelope beyond the first two oscillations
    ax_gamma.scatter(log_rho, Gamma, s=10, color=colors[idx],
                     label=f"B = {magnetic_fields[suffix]}")


# Load nonmagnetic BPS EOS data from CSV (format: rho,P,nB,Z,A,Gamma)
nonmagnetic_file = "/mnt/c/Users/karan/Dropbox/KARAN/2 Areas/Education/PhD/3 Research/pulsarmhd/data/bps.csv"
data_nonmag = pd.read_csv(nonmagnetic_file)

# Extract columns from nonmagnetic data
rho_nonmag = data_nonmag['rho'].values
P_nonmag = data_nonmag['P'].values
Gamma_nonmag = data_nonmag['Gamma'].values

# Compute logarithms for EOS plotting
log_rho_nonmag = np.log10(rho_nonmag)
log_P_nonmag = np.log10(P_nonmag)

# Optionally, filter out non-positive values (if necessary)
mask_nonmag = (rho_nonmag > 0.0) & (P_nonmag > 0.0)
log_rho_nonmag = log_rho_nonmag[mask_nonmag]
log_P_nonmag = log_P_nonmag[mask_nonmag]
Gamma_nonmag = Gamma_nonmag[mask_nonmag]

# Plot nonmagnetic EOS on the EOS figure with a distinct marker/color
ax_eos.scatter(log_rho_nonmag, log_P_nonmag, s=15, color='magenta', marker='x',
               label="Nonmagnetic BPS")

# Plot nonmagnetic Gamma on the Gamma figure
ax_gamma.scatter(log_rho_nonmag, Gamma_nonmag, s=15, color='magenta', marker='x',
                 label="Nonmagnetic BPS")

# Add legends to both plots
ax_eos.legend()
ax_gamma.legend()

# Save figures
fig_eos.savefig("magnetic_bps_eos_multi.png")
fig_gamma.savefig("magnetic_bps_gamma_multi.png")