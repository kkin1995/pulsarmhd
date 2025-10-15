import pandas as pd
import glob
import os

import matplotlib.pyplot as plt

M_sun = 1.989e33 # g

csv_files = glob.glob(os.path.join("/mnt/c/Users/karan/Dropbox/KARAN/2 Areas/Education/PhD/3 Research/gr_compact_objects_cpp/data", "*.csv"))

# Load the data
# filename = "electron_non_relativistic_rhoc_1.00pp05"
for file in csv_files:
    print(file)
    data = pd.read_csv(file)

    # Convert log values to actual values
    data['radius'] = 10 ** data['log_r[cm]']
    data['mass'] = 10 ** data['log_m[g]']
    data['pressure'] = 10 ** data['log_P[dyne/cm^2]']

    fig, ax1 = plt.subplots(figsize = (10, 6))

    ax1.scatter(data['radius'] * 1e-5, data['mass'] / M_sun, s=10, c='red', label='Mass')
    ax1.set_xlabel('Radius (r) [km]')
    ax1.set_ylabel('M / M_sun', color = 'red')
    ax1.tick_params(axis='y', labelcolor='red')

    ax2 = ax1.twinx()
    ax2.scatter(data['radius'] * 1e-5, data['pressure'], s=10, c='blue', label='Pressure')
    ax2.set_ylabel('Pressure (P) [dyn/cm^2]', color='blue')
    ax2.tick_params(axis='y', labelcolor='blue')

    plt.title('Radius vs Mass & Pressure')
    plt.legend()
    plt.grid(True)

    # # Add axes ticks
    # plt.xticks(np.arange(min(data['radius']), max(data['radius'])+1, 1.0))
    # plt.yticks(np.arange(min(min(data['mass']), min(data['pressure'])), max(max(data['mass']), max(data['pressure']))+1, 1.0))

    # Save the plot as a jpg file
    plt.savefig(os.path.join("plots", file, ".jpg"), format='jpg')
