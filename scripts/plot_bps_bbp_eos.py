import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Read the CSV files
bps_data = pd.read_csv('/mnt/c/Users/karan/Dropbox/KARAN/2 Areas/Education/PhD/3 Research/pulsarmhd/data/bps.csv')
bbp_data = pd.read_csv('/mnt/c/Users/karan/Dropbox/KARAN/2 Areas/Education/PhD/3 Research/pulsarmhd/data/bbp.csv')

# Find boundary region
bps_max = bps_data['rho'].max()
bps_min = bps_data['rho'].min()
bbp_max = bbp_data['rho'].max()
bbp_min = bbp_data['rho'].min()

# Calculate adiabatic index (Gamma) for each EOS
def calculate_gamma(density, pressure):
    # Convert to log space
    log_rho = np.log(density)
    log_p = np.log(pressure)
    
    # Calculate derivative using central differences
    gamma = np.gradient(log_p, log_rho)
    return gamma

# Calculate Gamma for both EOSs separately
bps_gamma = calculate_gamma(bps_data['rho'], bps_data['P'])
bbp_gamma = calculate_gamma(bbp_data['rho'], bbp_data['P'])

# Create figure with three subplots
fig = plt.figure(figsize=(15, 10))
gs = fig.add_gridspec(2, 2)

# Plot 1: Full range EOS
ax1 = fig.add_subplot(gs[0, :])
ax1.loglog(bps_data['rho'], bps_data['P'], 'b-', label='BPS EOS')
ax1.loglog(bbp_data['rho'], bbp_data['P'], 'r-', label='BBP EOS')
ax1.set_xlabel(r'$\rho$ (g/cm³)', fontsize=12)
ax1.set_ylabel(r'P (dyne/cm²)', fontsize=12)
ax1.set_title('Full Range EOS Comparison', fontsize=14)
ax1.grid(True, which="both", ls="-", alpha=0.2)
ax1.legend(fontsize=10)

# Print ranges
print("\nEOS Ranges:")
print(f"BPS: {bps_min:.2e} to {bps_max:.2e} g/cm³")
print(f"BBP: {bbp_min:.2e} to {bbp_max:.2e} g/cm³")

if bps_max < bbp_min:
    print("\nThere is a gap between the EOS")
    boundary_min = bps_max
    boundary_max = bbp_min
    plot_type = "gap"
else:
    print("\nThere is an overlap between the EOS")
    boundary_min = max(bps_min, bbp_min)
    boundary_max = min(bps_max, bbp_max)
    plot_type = "overlap"

# Plot 2: Zoomed boundary region
margin_factor = 2
plot_min = boundary_min / margin_factor
plot_max = boundary_max * margin_factor

ax2 = fig.add_subplot(gs[1, 0])
ax2.loglog(bps_data['rho'], bps_data['P'], 'b-', label='BPS EOS')
ax2.loglog(bbp_data['rho'], bbp_data['P'], 'r-', label='BBP EOS')
ax2.set_xlim(plot_min, plot_max)

# Find data in boundary region
bps_boundary = bps_data[(bps_data['rho'] >= plot_min) & (bps_data['rho'] <= plot_max)]
bbp_boundary = bbp_data[(bbp_data['rho'] >= plot_min) & (bbp_data['rho'] <= plot_max)]

if len(bps_boundary) > 0 and len(bbp_boundary) > 0:
    min_p = min(bps_boundary['P'].min(), bbp_boundary['P'].min()) / margin_factor
    max_p = max(bps_boundary['P'].max(), bbp_boundary['P'].max()) * margin_factor
    ax2.set_ylim(min_p, max_p)

ax2.set_xlabel(r'$\rho$ (g/cm³)', fontsize=12)
ax2.set_ylabel(r'P (dyne/cm²)', fontsize=12)
ax2.set_title(f'Boundary Region ({"Gap" if plot_type == "gap" else "Overlap"})', fontsize=14)
ax2.grid(True, which="both", ls="-", alpha=0.2)
ax2.legend(fontsize=10)

# Plot 3: Adiabatic index (avoiding the gap)
ax3 = fig.add_subplot(gs[1, 1])

# Plot BPS gamma up to slightly before the gap
bps_valid = bps_data['rho'] < bps_max * 0.99  # Stay slightly away from the edge
ax3.semilogx(bps_data.loc[bps_valid, 'rho'], bps_gamma[bps_valid], 'b-', label='BPS Γ')

# Plot BBP gamma starting slightly after the gap
bbp_valid = bbp_data['rho'] > bbp_min * 1.01  # Stay slightly away from the edge
ax3.semilogx(bbp_data.loc[bbp_valid, 'rho'], bbp_gamma[bbp_valid], 'r-', label='BBP Γ')

ax3.set_xlabel(r'$\rho$ (g/cm³)', fontsize=12)
ax3.set_ylabel(r'$\Gamma = d\log P/d\log \rho$', fontsize=12)
ax3.set_title('Adiabatic Index', fontsize=14)
ax3.grid(True, which="both", ls="-", alpha=0.2)
ax3.legend(fontsize=10)

# Set reasonable y-limits for gamma plot
ax3.set_ylim(0, 4)  # Adjust these values based on your data

# Print statistical summary of Gamma (excluding the gap region)
print("\nAdiabatic Index (Γ) Statistics:")
print("\nBPS EOS:")
print(f"Mean Γ: {np.mean(bps_gamma[bps_valid]):.3f}")
print(f"Min Γ: {np.min(bps_gamma[bps_valid]):.3f}")
print(f"Max Γ: {np.max(bps_gamma[bps_valid]):.3f}")
print(f"Std Γ: {np.std(bps_gamma[bps_valid]):.3f}")

print("\nBBP EOS:")
print(f"Mean Γ: {np.mean(bbp_gamma[bbp_valid]):.3f}")
print(f"Min Γ: {np.min(bbp_gamma[bbp_valid]):.3f}")
print(f"Max Γ: {np.max(bbp_gamma[bbp_valid]):.3f}")
print(f"Std Γ: {np.std(bbp_gamma[bbp_valid]):.3f}")

# Adjust layout
plt.tight_layout()

# Save the plot
plt.savefig('eos_comparison_with_gamma.png', dpi=300, bbox_inches='tight')