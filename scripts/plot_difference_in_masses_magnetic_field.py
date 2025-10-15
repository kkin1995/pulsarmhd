import pandas as pd
import glob
import math
import os
import matplotlib.pyplot as plt
import re
import numpy as np
from scipy.interpolate import interp1d

csv_pattern_hybrid_non_magnetic = '/mnt/c/Users/karan/Dropbox/KARAN/2 Areas/Education/PhD/3 Research/pulsarmhd/data/TOV-Solution-Non-Magnetic/tov_solution_magnetic_bps_bbp_polytrope_*.csv'
csv_pattern_hybrid_magnetic = '/mnt/c/Users/karan/Dropbox/KARAN/2 Areas/Education/PhD/3 Research/pulsarmhd/data/TOV-Solution-B-1/tov_solution_magnetic_bps_bbp_polytrope_*.csv'

file_list_hybrid_nonmagnetic = glob.glob(csv_pattern_hybrid_non_magnetic)
file_list_hybrid_magnetic = glob.glob(csv_pattern_hybrid_magnetic)

if not file_list_hybrid_nonmagnetic and not file_list_hybrid_magnetic:
    print("No CSV files found for either pattern. Check the file paths/patterns!")
    exit()

def get_last_valid_value(df, column):
    """
    Returns the last non-NaN value from the specified column in a DataFrame.
    """
    last_idx = df[column].last_valid_index()
    return df.loc[last_idx, column]

def parse_central_density_from_filename_hybrid(filename):
    """
    Parses the central density (in g/cm^3) from the hybrid EoS CSV filename.
    Expected pattern:
      .../tov_solution_magnetic_bps_bbp_polytrope_1.00e+16.csv
    Returns log10(rho_c) or None if not found.
    """
    base = os.path.basename(filename)
    match = re.search(r'polytrope_(.*)\.csv', base)
    if not match:
        return None
    rho_str = match.group(1)  # e.g. '1.00e+16'
    try:
        rho_c = float(rho_str)
        return math.log10(rho_c)
    except ValueError:
        return None

def extract_data(file_list, parse_func):
    """
    Reads each CSV in file_list, extracts the final (surface) mass and radius,
    plus parses log10(rho_c) from the filename using parse_func.
    Returns a list of tuples:
      [(mass_in_solar_masses, radius_in_km, log10_rho_c), ...]
    """
    data = []
    for filename in file_list:
        df = pd.read_csv(filename)
        # Check for required columns
        if 'log_m[g]' not in df.columns or 'log_r[cm]' not in df.columns:
            print(f"Skipping file {filename} due to missing columns.")
            continue

        # Get last valid log(m) and log(r)
        last_mass_log = get_last_valid_value(df, 'log_m[g]')
        last_radius_log = get_last_valid_value(df, 'log_r[cm]')

        # Convert from log scale to physical units
        mass_solar = 10**last_mass_log / 1.989e33  # from grams to M_sun
        radius_km  = 10**last_radius_log / 1e5     # from cm to km

        # Parse log10(rho_c)
        log_rho_c = parse_func(filename)
        if log_rho_c is None:
            print(f"Could not parse central density from filename: {filename}")
            continue

        data.append((mass_solar, radius_km, log_rho_c))
    return data

data_nonmagnetic = extract_data(file_list_hybrid_nonmagnetic, parse_central_density_from_filename_hybrid)
data_magnetic = extract_data(file_list_hybrid_magnetic, parse_central_density_from_filename_hybrid)

# Separate out each quantity for the hybrid EoS
masses_nonmagnetic = np.array([row[0] for row in data_nonmagnetic])
radii_nonmagnetic  = np.array([row[1] for row in data_nonmagnetic])
log_rho_cs_nonmagnetic = np.array([row[2] for row in data_nonmagnetic])

# Separate out each quantity for the neutron relativistic data
masses_magnetic = np.array([row[0] for row in data_magnetic])
radii_magnetic  = np.array([row[1] for row in data_magnetic])
log_rho_cs_magnetic = np.array([row[2] for row in data_magnetic])

relative_mass_difference = (masses_magnetic - masses_nonmagnetic) / masses_nonmagnetic

os.makedirs("plots", exist_ok=True)

plt.figure(figsize=(10, 6))

# Hybrid EoS
plt.scatter(radii_magnetic, relative_mass_difference, label='Hybrid EOS', color='blue', s=10)
plt.xlabel('Radius (km)')
plt.ylabel('Relative Mass Difference Between Magnetic BPS and Non-Magnetic BPS')
plt.title('Mass-Radius Relations')
plt.grid(True)
plt.legend()

plt.savefig(os.path.join("plots", "relative_mass_difference_radius_overplot.jpg"))
plt.clf()

plt.figure(figsize=(10, 6))

# Hybrid EoS
plt.scatter(log_rho_cs_magnetic, relative_mass_difference, label='Hybrid EOS', color='blue', s=10)
plt.xlabel('Radius (km)')
plt.ylabel('Relative Mass Difference Between Magnetic BPS and Non-Magnetic BPS')
plt.title('Mass - Central Density Relations')
plt.grid(True)
plt.legend()

plt.savefig(os.path.join("plots", "relative_mass_difference_density_overplot.jpg"))


# Create interpolation functions to compare at equal masses
# Sort the data by mass first
sorted_idx_nonmag = np.argsort(masses_nonmagnetic)
sorted_idx_mag = np.argsort(masses_magnetic)

masses_nonmag_sorted = masses_nonmagnetic[sorted_idx_nonmag]
radii_nonmag_sorted = radii_nonmagnetic[sorted_idx_nonmag]

masses_mag_sorted = masses_magnetic[sorted_idx_mag]
radii_mag_sorted = radii_magnetic[sorted_idx_mag]

# Create interpolation functions (use linear to avoid artifacts)
r_nonmag_interp = interp1d(masses_nonmag_sorted, radii_nonmag_sorted, kind='linear',
                          bounds_error=False, fill_value="extrapolate")
r_mag_interp = interp1d(masses_mag_sorted, radii_mag_sorted, kind='linear',
                       bounds_error=False, fill_value="extrapolate")

# Find common mass range
min_mass = max(np.min(masses_nonmagnetic), np.min(masses_magnetic))
max_mass = min(np.max(masses_nonmagnetic), np.max(masses_magnetic))

# Generate common mass points for comparison
common_masses = np.linspace(min_mass, max_mass, 100)

# Calculate radii at common mass points
radii_nonmag_interp = r_nonmag_interp(common_masses)
radii_mag_interp = r_mag_interp(common_masses)

# Calculate absolute and relative radius differences
delta_r = radii_mag_interp - radii_nonmag_interp
relative_r_diff = (radii_mag_interp - radii_nonmag_interp) / radii_nonmag_interp * 100  # as percentage

# Plot 2: Absolute radius difference vs mass
plt.figure(figsize=(10, 6))
plt.plot(common_masses, delta_r, 'g-', linewidth=2)
plt.axhline(y=0, color='k', linestyle='--', alpha=0.5)
plt.xlabel('Mass (Solar Masses)')
plt.ylabel('Radius Difference: Magnetic - Non-Magnetic (km)')
plt.title('Absolute Radius Difference vs. Mass')
plt.grid(True)
plt.savefig(os.path.join("plots", "absolute_radius_difference_vs_mass.jpg"))
plt.clf()

# Plot 3: Relative radius difference vs mass (percentage)
plt.figure(figsize=(10, 6))
plt.plot(common_masses, relative_r_diff, 'r-', linewidth=2)
plt.axhline(y=0, color='k', linestyle='--', alpha=0.5)
plt.xlabel('Mass (Solar Masses)')
plt.ylabel('Relative Radius Difference (%)')
plt.title('Percentage Radius Difference vs. Mass')
plt.grid(True)
plt.savefig(os.path.join("plots", "relative_radius_difference_percentage.jpg"))
plt.clf()

# Plot 4: Central density vs radius difference
plt.figure(figsize=(10, 6))
# Create interpolation functions for central density vs mass
# for both datasets to get central density at common mass points
log_rho_nonmag_sorted = log_rho_cs_nonmagnetic[sorted_idx_nonmag]
log_rho_mag_sorted = log_rho_cs_magnetic[sorted_idx_mag]

rho_nonmag_interp = interp1d(masses_nonmag_sorted, log_rho_nonmag_sorted, kind='linear',
                            bounds_error=False, fill_value="extrapolate")
rho_mag_interp = interp1d(masses_mag_sorted, log_rho_mag_sorted, kind='linear',
                          bounds_error=False, fill_value="extrapolate")

# Use the magnetic central densities for this plot
common_log_rho = rho_mag_interp(common_masses)

plt.scatter(common_log_rho, relative_r_diff, c=common_masses, s=10, cmap='viridis')
cbar = plt.colorbar()
cbar.set_label('Mass (Solar Masses)')
plt.xlabel('log10(Central Density) (g/cm³)')
plt.ylabel('Relative Radius Difference (%)')
plt.title('Radius Difference vs. Central Density')
plt.grid(True)
plt.savefig(os.path.join("plots", "radius_difference_vs_density.jpg"))

print("Analysis complete. Plots saved to 'plots' directory.")
