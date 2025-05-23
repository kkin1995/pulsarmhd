import pandas as pd
import matplotlib.pyplot as plt
import os

nonmagnetic_tov = pd.read_csv("/mnt/c/Users/karan/Dropbox/KARAN/2 Areas/Education/PhD/3 Research/pulsarmhd/data/TOV-Solution-Non-Magnetic/tov_solution_magnetic_bps_bbp_polytrope_1.00e+15.csv")
magnetic_tov = pd.read_csv("/mnt/c/Users/karan/Dropbox/KARAN/2 Areas/Education/PhD/3 Research/pulsarmhd/data/TOV-Solution-B-1/tov_solution_magnetic_bps_bbp_polytrope_1.00e+15.csv")

# log_r[cm],log_m[g],log_P[dyne/cm^2]

log_r_nonmagnetic = nonmagnetic_tov['log_r[cm]']
log_m_nonmagnetic = nonmagnetic_tov['log_m[g]']
log_P_nonmagnetic = nonmagnetic_tov['log_P[dyne/cm^2]']

log_r_magnetic = magnetic_tov['log_r[cm]']
log_m_magnetic = magnetic_tov['log_m[g]']
log_P_magnetic = magnetic_tov['log_P[dyne/cm^2]']

nonmagnetic_mass = 10**log_m_nonmagnetic
nonmagnetic_mass_solar = 10**log_m_nonmagnetic / 1.989e33  # from grams to M_sun
nonmagnetic_radius_km  = 10**log_r_nonmagnetic / 1e5     # from cm to km

magnetic_mass = 10**log_m_magnetic
magnetic_mass_solar = 10**log_m_magnetic / 1.989e33  # from grams to M_sun
magnetic_radius_km  = 10**log_r_magnetic / 1e5     # from cm to km

relative_mass_difference = ((nonmagnetic_mass - magnetic_mass) / nonmagnetic_mass) * 100

plt.figure(figsize=(10, 6))

# Hybrid EoS
plt.scatter(nonmagnetic_radius_km, relative_mass_difference, label='Hybrid EOS', color='blue', s=10)
plt.xlabel('Radius (km)')
plt.ylabel('Percentage Relative Mass Difference Between Magnetic BPS and Non-Magnetic BPS')
plt.title('Mass-Radius Relations')
plt.grid(True)
plt.legend()

plt.savefig(os.path.join("plots", "relative_mass_difference_radius_1.00E15.jpg"))
plt.clf()
