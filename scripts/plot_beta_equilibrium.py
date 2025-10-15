import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Read the data
df = pd.read_csv('/mnt/c/Users/karan/Dropbox/KARAN/2 Areas/Education/PhD/3 Research/pulsarmhd/beta_equilibrium_function.csv')

# Create the plot
plt.figure(figsize=(10, 6))
plt.plot(df['x_n'], df['f_value'], 'b-', label='β-equilibrium function')

# Add horizontal line at y=0
plt.axhline(y=0, color='r', linestyle='--', label='f(x) = 0')

# Add labels and title
plt.xlabel('x_n')
plt.ylabel('f(x_n)')
plt.title('Beta Equilibrium Function vs x_n')
plt.grid(True)
plt.legend()

# Scientific notation for y-axis
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

# Save the plot
plt.savefig('beta_equilibrium_plot.png', dpi=300, bbox_inches='tight')
