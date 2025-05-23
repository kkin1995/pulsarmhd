import pandas as pd
import matplotlib.pyplot as plt

# Load data from CSV
data_file = "/mnt/c/Users/karan/Dropbox/KARAN/2 Areas/Education/PhD/3 Research/pulsarmhd/data/tov_solution_magnetic_bps_bbp_polytrope_7.50e+17.csv"
df = pd.read_csv(data_file)
df.dropna(inplace=True)

# Extract logarithmic data for radius, mass, and pressure
log_r = df["log_r[cm]"].values       # Logarithm of radius in cm
log_m = df["log_m[g]"].values         # Logarithm of mass in g
log_P = df["log_P[dyne/cm^2]"].values  # Logarithm of pressure in dyne/cm^2

# Convert from logarithmic to linear scale
r_cm = 10 ** log_r      # Radius in cm
m_g = 10 ** log_m       # Mass in g
p = 10 ** log_P         # Pressure in dyne/cm^2

# Unit conversions:
r_km = r_cm / 1e5       # Convert radius from cm to km
M_sun = 1.989e33        # Solar mass in g
m_solar = m_g / M_sun   # Mass in solar masses

# Create figure and twin axes for dual y-axis plotting
fig, ax1 = plt.subplots(figsize=(10, 5))
ax2 = ax1.twinx()

# Plot mass on the primary (left) y-axis
ax1.plot(r_km, m_solar, 'r-', label='Mass Profile')
ax1.set_xlabel("Radius [km]")
ax1.set_ylabel("Mass [M_sun]", color='red')
ax1.tick_params(axis='y', colors='red')

# Plot pressure on the secondary (right) y-axis
ax2.plot(r_km, p, 'b-', label='Pressure Profile')
ax2.set_ylabel("Pressure [dyne/cm^2]", color='blue')
ax2.tick_params(axis='y', colors='blue')

# Add title and grid
plt.title("Mass and Pressure Profiles")
ax1.grid(True)

# Optionally, create a combined legend
lines_1, labels_1 = ax1.get_legend_handles_labels()
lines_2, labels_2 = ax2.get_legend_handles_labels()
ax1.legend(lines_1 + lines_2, labels_1 + labels_2, loc='best')

plt.tight_layout()
plt.savefig("plots/tov_m_r_plot_magnetic_bps_bbp_neutron_polytrope_7.50e+17.png")
