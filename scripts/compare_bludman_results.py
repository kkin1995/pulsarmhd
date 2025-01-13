import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import PchipInterpolator
import os
import glob

G = 6.67e-8 # CGS
c = 3e10 # CGS

def convert_bludman_data_to_mass_and_radius(path_to_bludman_results: str, n: float, K: float) -> tuple[np.ndarray, np.ndarray]:
    """
    Convert Bludman's dimensionless quantities to physical mass and radius.
    
    Args:
        path_to_bludman_results: Path to CSV file containing Bludman's results
        n: Polytropic index
        K: Polytropic constant in CGS units
    
    Returns:
        Tuple of (masses in solar masses, radii in km)
    """
    df = pd.read_csv(path_to_bludman_results)

    sigma = df.loc[:, "sigma"].values
    zeta_1 = df.loc[:, "zeta_1"].values
    v_zeta_1 = df.loc[:, "v(zeta_1)"].values
    M_tilde = df.loc[:, "M_tilde"].values

    mask = sigma != 0
    sigma = sigma[mask]
    zeta_1 = zeta_1[mask]
    v_zeta_1 = v_zeta_1[mask]
    M_tilde = M_tilde[mask]

    A = (((4 * np.pi * G) / (c ** 2)) * (1 / (n + 1)) * ((c **2 / K) ** n) * sigma ** (n - 1)) ** (1/2)

    # Radius
    R = (1 / A) * zeta_1

    # Mass
    M = ( (1 / (4 * np.pi)) * ( ((n + 1) * c**2) / G )**3 * (K / c**2)**n)**(1/2) * M_tilde

    mass_in_solar_mass = M / 1.989e33
    radius_in_kms = R / 1e5

    return mass_in_solar_mass, radius_in_kms

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

def plot_mass_radius(values, label, color):
    masses, radii = zip(*values)
    masses = [10**mass / 1.989e33 for mass in masses]  # Convert mass from grams to solar masses
    radii = [10**radius / 1e5 for radius in radii]  # Convert radius from cm to km
    plt.scatter(radii, masses, label=label, color=color, s = 10)

if __name__ == "__main__":
    neutron_relativistic_dfs = [pd.read_csv(file) for file in glob.glob('/mnt/c/Users/karan/Dropbox/KARAN/2 Areas/Education/PhD/3 Research/pulsarmhd/data/neutron_relativistic*.csv')]
    neutron_relativistic_values = extract_last_values(neutron_relativistic_dfs)
    simulated_mass, simulated_radii = zip(*neutron_relativistic_values)
    masses = [10**mass / 1.989e33 for mass in simulated_mass]  # Convert mass from grams to solar masses
    radii = [10**radius / 1e5 for radius in simulated_radii]  # Convert radius from cm to km

    # Relativistic Neutron Gas
    n = 3
    K = 1.2293e15
    path_to_bludman_results = "/mnt/c/Users/karan/Dropbox/KARAN/2 Areas/Education/PhD/3 Research/pulsarmhd/data/bludman_results/n=3-and-sigma_CR=0.csv"
    mass_in_solar_mass, radius_in_kms = convert_bludman_data_to_mass_and_radius(path_to_bludman_results, n, K)

    plt.figure(figsize=(10, 6))
    plt.xlabel('Radius (km)')
    plt.ylabel('Mass (Solar Masses)')
    plot_mass_radius(neutron_relativistic_values, "Karan Simulations", color="blue")
    plt.scatter(radius_in_kms, mass_in_solar_mass, label="Bludman Paper", color="red", s=10)
    plt.grid(True)
    plt.legend()
    plt.title('Relativistic Neutron Mass-Radius Relations')
    plt.savefig(os.path.join("/mnt/c/Users/karan/Dropbox/KARAN/2 Areas/Education/PhD/3 Research/pulsarmhd/data/bludman_results", "n=3-and-sigma_CR=0.jpg"))