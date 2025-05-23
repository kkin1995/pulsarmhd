import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

# Define file paths
bbp_file_path = "/mnt/c/Users/karan/Dropbox/KARAN/2 Areas/Education/PhD/3 Research/pulsarmhd/data/bbp.csv"

# Define magnetic field values and corresponding file names
magnetic_fields = ["1", "01", "001", "10", "100"]
bps_file_template = "/mnt/c/Users/karan/Dropbox/KARAN/2 Areas/Education/PhD/3 Research/pulsarmhd/magnetic_bps_eos_B_{}.csv"

# Load BBP data
bbp_data = pd.read_csv(bbp_file_path)

# Set up plot for EoS curves
plt.figure(figsize=(8, 6))

# Plot BBP EoS
plt.plot(
    np.log10(bbp_data["rho"]), np.log10(bbp_data["P"]),
    label="BBP EoS", linestyle="-", marker="o", markersize=3
)

# Loop through different Magnetic BPS files and plot them
for B in magnetic_fields:
    bps_file_path = bps_file_template.format(B)
    
    try:
        # Load Magnetic BPS data
        bps_data = pd.read_csv(bps_file_path)

        # Plot Magnetic BPS (directly using log data)
        plt.plot(
            bps_data["log_rho"], bps_data["log_P"],
            label=f"Magnetic BPS EoS (B = {B} B_c)", linestyle="--", marker="s", markersize=3
        )
    except FileNotFoundError:
        print(f"Warning: File not found for B = {B} B_c. Skipping.")

# Labels and title
plt.xlabel(r"Log Density $\log_{10}(\rho)$ (g/cm$^3$)", fontsize=14)
plt.ylabel(r"Log Pressure $\log_{10}(P)$ (dyne/cm$^2$)", fontsize=14)
plt.title("Equation of State: BBP vs. Magnetic BPS (Various B Fields)", fontsize=16)
plt.legend(fontsize=10)

# Grid and academic styling
plt.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.7)
plt.tick_params(axis='both', which='major', labelsize=12)

# Save plot
plt.savefig("bbp_magnetic_bps_all.png")

# Clear the figure
plt.clf()

# Set up plot for pressure differences
plt.figure(figsize=(8, 6))

# Compute pressure difference for each BPS file
for B in magnetic_fields:
    bps_file_path = bps_file_template.format(B)
    
    try:
        # Load Magnetic BPS data
        bps_data = pd.read_csv(bps_file_path)

        # Create interpolation function for Magnetic BPS pressure
        bps_interp = interp1d(bps_data["log_rho"], bps_data["log_P"], kind="linear", fill_value="extrapolate")

        # Compute pressure differences at BBP densities
        bbp_log_rho = np.log10(bbp_data["rho"])
        bbp_log_P = np.log10(bbp_data["P"])
        bps_log_P_interp = bps_interp(bbp_log_rho)  # Interpolated BPS pressure at BBP densities

        pressure_diff = bbp_log_P - bps_log_P_interp  # Difference P_BBP - P_BPS

        # Plot the pressure difference
        plt.plot(bbp_log_rho, pressure_diff, label=f"B = {B} B_c", linestyle="-", marker="o", markersize=2)

    except FileNotFoundError:
        print(f"Warning: File not found for B = {B} B_c. Skipping.")

# Indicate zero crossing
plt.axhline(0, color='k', linestyle="--", linewidth=1)

# Labels and title
plt.xlabel(r"Log Density $\log_{10}(\rho)$ (g/cm$^3$)", fontsize=14)
plt.ylabel(r"Pressure Difference $\log_{10}(P_{\rm BBP}) - \log_{10}(P_{\rm BPS})$", fontsize=14)
plt.title("Pressure Difference: BBP vs. Magnetic BPS (Various B Fields)", fontsize=16)
plt.legend(fontsize=10)

# Grid and academic styling
plt.grid(True, linestyle="--", linewidth=0.5, alpha=0.7)
plt.tick_params(axis='both', which='major', labelsize=12)

# Save plot
plt.savefig("bbp_magnetic_bps_pressure_diff_all.png")
