import pandas as pd
import glob
import os
import matplotlib.pyplot as plt
import math
import re

# ------------------------------------------------------------------------
# 1. Gather CSV files for the hybrid EoS
# ------------------------------------------------------------------------
csv_pattern_hybrid = '/mnt/c/Users/karan/Dropbox/KARAN/2 Areas/Education/PhD/3 Research/pulsarmhd/data/tov_solution_magnetic_bps_bbp_polytrope_*.csv'
file_list_hybrid = glob.glob(csv_pattern_hybrid)

# Gather CSV files for the neutron relativistic data
csv_pattern_rel = '/mnt/c/Users/karan/Dropbox/KARAN/2 Areas/Education/PhD/3 Research/pulsarmhd/data/neutron_relativistic_rhoc_*.csv'
file_list_rel = glob.glob(csv_pattern_rel)

if not file_list_hybrid and not file_list_rel:
    print("No CSV files found for either pattern. Check the file paths/patterns!")
    exit()

# ------------------------------------------------------------------------
# 2. Helper functions
# ------------------------------------------------------------------------
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

def parse_central_density_from_filename_rel(filename):
    """
    Parses the central density (in g/cm^3) from the neutron relativistic CSV filename.
    Expected pattern:
      .../neutron_relativistic_rhoc_5.00pp18.csv
    which corresponds to 5.00e+18 in standard notation.
    Returns log10(rho_c) or None if not found.
    """
    base = os.path.basename(filename)
    match = re.search(r'neutron_relativistic_rhoc_(.*)\.csv', base)
    if not match:
        return None
    rho_str = match.group(1)  # e.g. '5.00pp18'
    # Convert '5.00pp18' -> '5.00e+18'
    rho_str = rho_str.replace('pp', 'e+')
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

# ------------------------------------------------------------------------
# 3. Extract data for both sets
# ------------------------------------------------------------------------
data_hybrid = extract_data(file_list_hybrid, parse_central_density_from_filename_hybrid)
data_rel = extract_data(file_list_rel, parse_central_density_from_filename_rel)

# Separate out each quantity for the hybrid EoS
masses_hybrid = [row[0] for row in data_hybrid]
radii_hybrid  = [row[1] for row in data_hybrid]
log_rho_cs_hybrid = [row[2] for row in data_hybrid]

# Separate out each quantity for the neutron relativistic data
masses_rel = [row[0] for row in data_rel]
radii_rel  = [row[1] for row in data_rel]
log_rho_cs_rel = [row[2] for row in data_rel]

if not (masses_hybrid or masses_rel):
    print("No valid data found in the CSV files.")
    exit()

os.makedirs("plots", exist_ok=True)

# ------------------------------------------------------------------------
# 4. Plot 1: Mass vs. Radius
# ------------------------------------------------------------------------
plt.figure(figsize=(10, 6))

# Hybrid EoS
if masses_hybrid:
    plt.scatter(radii_hybrid, masses_hybrid, label='Hybrid EOS', color='blue', s=10)

# Neutron Relativistic
if masses_rel:
    plt.scatter(radii_rel, masses_rel, label='Neutron Relativistic', color='green', s=10)

plt.xlabel('Radius (km)')
plt.ylabel('Mass (Solar Masses)')
plt.title('Mass-Radius Relations')
plt.grid(True)
plt.legend()

plt.savefig(os.path.join("plots", "mass_radius_overplot.jpg"))
plt.clf()

# ------------------------------------------------------------------------
# 5. Plot 2: Mass vs. log10(rho_c)
# ------------------------------------------------------------------------
plt.figure(figsize=(10, 6))

# Hybrid EoS
if masses_hybrid:
    plt.scatter(log_rho_cs_hybrid, masses_hybrid, label='Hybrid EOS', color='red', s=10)

# Neutron Relativistic
if masses_rel:
    plt.scatter(log_rho_cs_rel, masses_rel, label='Neutron Relativistic', color='purple', s=10)

plt.xlabel('log10(rho_c) [g/cm^3]')
plt.ylabel('Mass (Solar Masses)')
plt.title('Mass vs. log10(Central Density)')
plt.grid(True)
plt.legend()

plt.savefig(os.path.join("plots", "mass_vs_logrho_c_overplot.jpg"))
plt.clf()
