import pandas as pd
import glob
import os
import matplotlib.pyplot as plt

electron_non_relativistic_dfs = [pd.read_csv(file) for file in glob.glob('/mnt/c/Users/karan/Dropbox/KARAN/2 Areas/Education/PhD/3 Research/gr_compact_objects_cpp/data/electron_non_relativistic*.csv')]
electron_relativistic_dfs = [pd.read_csv(file) for file in glob.glob('/mnt/c/Users/karan/Dropbox/KARAN/2 Areas/Education/PhD/3 Research/gr_compact_objects_cpp/data/electron_relativistic*.csv')]
neutron_non_relativistic_dfs = [pd.read_csv(file) for file in glob.glob('/mnt/c/Users/karan/Dropbox/KARAN/2 Areas/Education/PhD/3 Research/gr_compact_objects_cpp/data/neutron_non_relativistic*.csv')]
neutron_relativistic_dfs = [pd.read_csv(file) for file in glob.glob('/mnt/c/Users/karan/Dropbox/KARAN/2 Areas/Education/PhD/3 Research/gr_compact_objects_cpp/data/neutron_relativistic*.csv')]

def get_last_valid_value(df, column):
    if pd.isna(df[column].iloc[-1]):
        return df[column].iloc[-2]
    return df[column].iloc[-1]

def extract_last_values(dfs):
    last_values = []
    for df in dfs:
        last_mass = get_last_valid_value(df, 'log_m[g]')
        last_radius = get_last_valid_value(df, 'log_r[cm]')
        last_values.append((last_mass, last_radius))
    return last_values

electron_non_relativistic_values = extract_last_values(electron_non_relativistic_dfs)
electron_relativistic_values = extract_last_values(electron_relativistic_dfs)
neutron_non_relativistic_values = extract_last_values(neutron_non_relativistic_dfs)
neutron_relativistic_values = extract_last_values(neutron_relativistic_dfs)

def plot_mass_radius(values, label, color):
    masses, radii = zip(*values)
    masses = [10**mass / 1.989e33 for mass in masses]  # Convert mass from grams to solar masses
    radii = [10**radius / 1e5 for radius in radii]  # Convert radius from cm to km
    plt.scatter(radii, masses, label=label, color=color, s = 10)

plt.figure(figsize=(10, 6))

plt.xlabel('Radius (km)')
plt.ylabel('Mass (Solar Masses)')
plot_mass_radius(electron_non_relativistic_values, 'Electron Non-Relativistic', 'blue')
plot_mass_radius(electron_relativistic_values, 'Electron Relativistic', 'red')
plt.ylim(-1, 2)
plt.legend()
plt.title('Mass-Radius Relations')
plt.grid(True)
plt.savefig(os.path.join("plots", "electron_gas_mass_radius.jpg"))
plt.clf()

plt.figure(figsize=(10, 6))

plt.xlabel('Radius (km)')
plt.ylabel('Mass (Solar Masses)')
plot_mass_radius(neutron_non_relativistic_values, 'Neutron Non-Relativistic', 'green')
plot_mass_radius(neutron_relativistic_values, 'Neutron Relativistic', 'purple')
plt.ylim(-1, 7)
plt.legend()
plt.title('Mass-Radius Relations')
plt.grid(True)
plt.savefig(os.path.join("plots", "neutron_gas_mass_radius.jpg"))
plt.clf()