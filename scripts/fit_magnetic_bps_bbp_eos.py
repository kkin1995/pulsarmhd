import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import math

def fritsch_carlson_monotonic_slopes(x1, x2, y1, y2, m1, m2):
    """
    Apply Fritsch-Carlson monotonicity constraints to cubic Hermite spline derivatives.
    
    Parameters:
    ----------
    x1, y1 : float
        First point coordinates (log_rho, log_P for EOS)
    x2, y2 : float
        Second point coordinates (log_rho, log_P for EOS)
    m1, m2 : float
        Initial derivatives at the first and second points
        
    Returns:
    -------
    tuple
        Adjusted derivatives (m1_new, m2_new) ensuring monotonicity
    """
    # Step 1: Calculate secant slope
    dx = x2 - x1
    dy = y2 - y1
    
    if abs(dx) < 1e-10:  # Vertical line
        return 0.0, 0.0
    
    delta = dy / dx  # Secant slope
    
    # Step 2: Check if delta is zero or near-zero (flat segment)
    if abs(delta) < 1e-10:
        return 0.0, 0.0
    
    # Step 3: Calculate non-dimensional ratios
    alpha = m1 / delta
    beta = m2 / delta
    
    # Step 4: Check if slopes have opposite signs from secant
    if alpha < 0 or beta < 0:
        # Different signs - not monotonic, set both to zero
        return 0.0, 0.0
    
    # Step 5: Check the three monotonicity conditions
    # Condition 1: if alpha*(2*alpha + beta - 3)^2/(3*(alpha + beta - 2)) > 0
    # Condition 2: if alpha + 2*beta - 3 <= 0
    # Condition 3: if 2*alpha + beta - 3 <= 0
    
    # Calculate terms for conditions
    sum_ab = alpha + beta
    cond2 = alpha + 2*beta - 3
    cond3 = 2*alpha + beta - 3
    
    # Check if at least one condition is satisfied
    if sum_ab == 2:  # Special case: linear
        # Already monotonic if endpoints are monotonic
        return m1, m2
    
    cond1_satisfied = False
    if sum_ab > 2:  # Parabola opens upward
        # Complex condition 1 check
        term1 = 2*alpha + beta - 3
        term2 = term1*term1/(3*(sum_ab - 2))
        cond1_satisfied = alpha*term2 > 0
    
    # If any condition is satisfied, original slopes are fine
    if cond1_satisfied or cond2 <= 0 or cond3 <= 0:
        return m1, m2
    
    # Step 6: No condition satisfied - apply scaling to ensure monotonicity
    # Simple approach: ensure alpha + beta <= 3 (a known sufficient condition)
    if sum_ab > 3:
        # Scale both derivatives
        tau = 3.0 / sum_ab
        return tau * m1, tau * m2
    
    # Alternative: apply circle constraint as a general safety measure
    alpha_beta_squared = alpha*alpha + beta*beta
    if alpha_beta_squared > 9:
        tau = 3.0 / math.sqrt(alpha_beta_squared)
        return tau * m1, tau * m2
    
    # If we reach here, the slopes should be fine, but use smaller values for safety
    # Reduce both by a safety factor as a last resort
    safety_factor = 0.9  # Slightly reduce derivatives
    return safety_factor * m1, safety_factor * m2

def linear_interpolation(x, x1, x2, y1, y2):
    t = (x - x1) / (x2 - x1)
    return y1 * (1-t) + y2 * t

# READING DATA
bbp_data = pd.read_csv("/mnt/c/Users/karan/Dropbox/KARAN/2 Areas/Education/PhD/3 Research/pulsarmhd/data/bbp.csv")
magnetic_bps_data = pd.read_csv("/mnt/c/Users/karan/Dropbox/KARAN/2 Areas/Education/PhD/3 Research/pulsarmhd/data/magnetic_bps_eos_B_001.csv")
bps_data = pd.read_csv("/mnt/c/Users/karan/Dropbox/KARAN/2 Areas/Education/PhD/3 Research/pulsarmhd/data/bps.csv")

# rho,P,nb,A,Z,Gamma
log_rho_bbp = np.log10(bbp_data["rho"].values)
log_P_bbp = np.log10(bbp_data["P"].values)

# log_n,log_rho,log_P
log_rho_bps = magnetic_bps_data["log_rho"].values
log_P_bps = magnetic_bps_data["log_P"].values

# # rho,P,nB,Z,A,Gamma - Non-Magnetic BPS EOS
# rho = bps_data["rho"].values
# P = bps_data["P"].values
# log_rho_bps = np.log10(rho)
# log_P_bps = np.log10(P)

# REMOVING NAN VALUES
# Create masks for valid values
valid_bps_mask = ~np.isnan(log_rho_bps) & ~np.isnan(log_P_bps) & ~np.isinf(log_rho_bps) & ~np.isinf(log_P_bps)
valid_bbp_mask = ~np.isnan(log_rho_bbp) & ~np.isnan(log_P_bbp) & ~np.isinf(log_rho_bbp) & ~np.isinf(log_P_bbp)

# Apply masks to both arrays simultaneously
log_rho_bps = log_rho_bps[valid_bps_mask]
log_P_bps = log_P_bps[valid_bps_mask]
log_rho_bbp = log_rho_bbp[valid_bbp_mask]
log_P_bbp = log_P_bbp[valid_bbp_mask]

# FITTING CUBIC SPLINES TO FIND POINT OF CLOSEST APPROACH OF BOTH EOS's
log_rho_min = np.max([np.min(log_rho_bps), np.min(log_rho_bbp)])
log_rho_max = np.min([np.max(log_rho_bps), np.max(log_rho_bbp)])
log_rho_common = np.linspace(log_rho_min, log_rho_max, 1000)
unique_log_rho_bps, indices = np.unique(log_rho_bps, return_index=True)
unique_log_P_bps = log_P_bps[indices]
print(f"BPS density range: {np.min(unique_log_rho_bps)} to {np.max(unique_log_rho_bps)}")
print(f"BBP density range: {np.min(log_rho_bbp)} to {np.max(log_rho_bbp)}")
f_bps = interp1d(unique_log_rho_bps, unique_log_P_bps, kind='cubic', bounds_error=False, fill_value="extrapolate")
f_bbp = interp1d(log_rho_bbp, log_P_bbp, kind='cubic', bounds_error=False, fill_value="extrapolate")
log_P_bps_interpolated = f_bps(log_rho_common)
log_P_bbp_interpolated = f_bbp(log_rho_common)
distance_bps_bbp = np.abs(log_P_bps_interpolated - log_P_bbp_interpolated)
idx_min_distance = np.argmin(distance_bps_bbp)
log_rho_min_distance = log_rho_common[idx_min_distance]
print(r'$\log_{10}(\rho)$ at minimum distance: ', log_rho_min_distance)

log_P_bps_transition = f_bps(log_rho_min_distance)
log_P_bbp_transition = f_bbp(log_rho_min_distance)
print(f"log(P) from BPS at transition: {log_P_bps_transition}")
print(f"log(P) from BBP at transition: {log_P_bbp_transition}")
print(f"Pressure difference: {10 ** log_P_bps_transition - 10 ** log_P_bbp_transition}")

def phi0(t):
    return 1 - (3 * (t ** 2)) + (2 * (t ** 3))

def phi1(t):
    return (3 * (t ** 2)) - (2 * (t ** 3))

def psi0(t):
    return (t ** 3) - (2 * (t ** 2)) + t

def psi1(t):
    return (t ** 3) - (t ** 2)

delta_1 = 0.2
delta_2 = 0.35

# Non-Magnetic BPS
# delta_1 = 0.0001
# delta_2 = 0.0001

x_1 = log_rho_min_distance - delta_1
x_2 = log_rho_min_distance + delta_2
h = x_2 - x_1 # Transition Width
f_1 = f_bps(x_1) # Value of Magnetic BPS Pressure at x_1
f_2 = f_bbp(x_2) # Value of BBP Pressure at x_2

# Define small delta for numerical derivative
delta = 0.001

# Calculate adiabatic index for BPS at x₁ (transition_start)
# Use forward difference
f_1_prime = (f_bps(x_1 + delta) - f_bps(x_1)) / delta

# Calculate adiabatic index for BBP at x₂ (transition_end)
# Since we're near the end of data range, use backward difference to be safe
f_2_prime = (f_bbp(x_2) - f_bbp(x_2 - delta)) / delta

def cubic_spline(x, f_1, f_2, f_1_prime, f_2_prime):
    t = (x - x_1) / (x_2 - x_1)
    first_term = f_1 * phi0(t)
    second_term = f_2 * phi1(t)
    third_term = h * f_1_prime * psi0(t)
    fourth_term = h * f_2_prime * psi1(t)

    return first_term + second_term + third_term + fourth_term

# Create dense grid in transition region
transition_x = np.linspace(x_1, x_2, 100)
transition_y = np.array([cubic_spline(x, f_1, f_2, f_1_prime, f_2_prime) for x in transition_x])

# Select BPS data below transition
bps_mask = unique_log_rho_bps < x_1
bps_hybrid_rho = unique_log_rho_bps[bps_mask]
bps_hybrid_P = unique_log_P_bps[bps_mask]

# Select BBP data above transition
bbp_mask = log_rho_bbp > x_2
bbp_hybrid_rho = log_rho_bbp[bbp_mask]
bbp_hybrid_P = log_P_bbp[bbp_mask]

# Combine all regions
hybrid_log_rho = np.concatenate([bps_hybrid_rho, transition_x, bbp_hybrid_rho])
hybrid_log_P = np.concatenate([bps_hybrid_P, transition_y, bbp_hybrid_P])

# Sort by density to ensure proper ordering
sort_idx = np.argsort(hybrid_log_rho)
hybrid_log_rho = hybrid_log_rho[sort_idx]
hybrid_log_P = hybrid_log_P[sort_idx]

# ---

rho_neutron_polytrope = np.linspace(1e12, 1e15, 100)
neutron_polytrope_k = 5.3802e9 # CGS
neutron_polytrope_gamma = 5.0 / 3.0
neutron_polytrope_P = neutron_polytrope_k * (rho_neutron_polytrope ** neutron_polytrope_gamma)

log_neutron_polytrope_rho = np.log10(rho_neutron_polytrope)[1:]
log_neutron_polytrope_P = np.log10(neutron_polytrope_P)[1:]

f_hybrid = interp1d(hybrid_log_rho, hybrid_log_P, kind='cubic', bounds_error=False, fill_value='extrapolate')

# 2. Find transition point between hybrid EOS and polytrope
log_rho_hybrid_max = np.max(hybrid_log_rho)
log_rho_poly_min = np.min(log_neutron_polytrope_rho)

if log_rho_hybrid_max >= log_rho_poly_min:
    # Overlap exists - find minimum distance point
    log_rho_overlap = np.linspace(log_rho_poly_min, log_rho_hybrid_max, 1000)
    log_P_hybrid_overlap = f_hybrid(log_rho_overlap)
    log_P_poly_overlap = np.log10(neutron_polytrope_k * (10**log_rho_overlap)**neutron_polytrope_gamma)
    
    distance_hybrid_poly = np.abs(log_P_hybrid_overlap - log_P_poly_overlap)
    idx_min_distance = np.argmin(distance_hybrid_poly)
    log_rho_transition = log_rho_overlap[idx_min_distance]
    
    print(f"Found overlap between hybrid EOS and polytrope")
    print(f"Transition point at log(rho) = {log_rho_transition:.4f}")
else:
    # No overlap - use midpoint between datasets
    log_rho_transition = (log_rho_hybrid_max + log_rho_poly_min) / 2
    print(f"No overlap between hybrid EOS and polytrope")
    print(f"Using midpoint for transition: log(rho) = {log_rho_transition:.4f}")

# 3. Define transition region
poly_delta_1 = 0.1  # Width before transition
poly_delta_2 = 0.1  # Width after transition
mid_rho = (log_rho_hybrid_max + log_rho_poly_min) / 2
poly_x_1 = log_rho_hybrid_max - poly_delta_1
poly_x_2 = log_rho_poly_min + poly_delta_2
poly_h = poly_x_2 - poly_x_1  # Transition width

# 4. Calculate function values at transition boundaries
poly_f_1 = f_hybrid(poly_x_1)  # Hybrid pressure at transition start
# Polytrope pressure at transition end
poly_f_2 = np.log10(neutron_polytrope_k) + poly_x_2 * neutron_polytrope_gamma

# 5. Calculate derivatives at transition boundaries
delta = 0.001  # Small step for numerical derivative
poly_f_1_prime = (f_hybrid(poly_x_1 + delta) - f_hybrid(poly_x_1)) / delta
poly_f_2_prime = neutron_polytrope_gamma  # Analytical derivative for polytrope

# print(f"Before adjustment - Derivatives: {poly_f_1_prime:.4f}, {poly_f_2_prime:.4f}")

# poly_f_1_prime, poly_f_2_prime = fritsch_carlson_monotonic_slopes(
#     poly_x_1, poly_x_2, poly_f_1, poly_f_2, poly_f_1_prime, poly_f_2_prime)

# print(f"After adjustment - Derivatives: {poly_f_1_prime:.4f}, {poly_f_2_prime:.4f}")

# Add this right before creating the transition region
print("\nDEBUGGING BBP-POLYTROPE TRANSITION:")
print(f"Transition point 1: log(rho)={poly_x_1:.4f}, log(P)={poly_f_1:.4f}")
print(f"Transition point 2: log(rho)={poly_x_2:.4f}, log(P)={poly_f_2:.4f}")
print(f"Pressure difference: {poly_f_2 - poly_f_1:.4f} (in log space)")
print(f"Derivative at point 1: {poly_f_1_prime:.4f}")
print(f"Derivative at point 2: {poly_f_2_prime:.4f}")

# 6. Create transition region using cubic spline
poly_transition_x = np.linspace(poly_x_1, poly_x_2, 100)
# poly_transition_y = np.array([cubic_spline(x, poly_f_1, poly_f_2, 
#                                            poly_f_1_prime, poly_f_2_prime) 
#                               for x in poly_transition_x])

poly_transition_y = np.array([linear_interpolation(x, poly_x_1, poly_x_2, 
                                                  poly_f_1, poly_f_2) 
                              for x in poly_transition_x])

# 7. Select data for final unified EOS
hybrid_mask = hybrid_log_rho < poly_x_1
final_hybrid_rho = hybrid_log_rho[hybrid_mask]
final_hybrid_P = hybrid_log_P[hybrid_mask]

poly_mask = log_neutron_polytrope_rho > poly_x_2
final_poly_rho = log_neutron_polytrope_rho[poly_mask]
final_poly_P = log_neutron_polytrope_P[poly_mask]

# 8. Combine all regions into unified EOS
unified_log_rho = np.concatenate([final_hybrid_rho, poly_transition_x, final_poly_rho])
unified_log_P = np.concatenate([final_hybrid_P, poly_transition_y, final_poly_P])

# 9. Sort by density to ensure proper ordering
sort_idx = np.argsort(unified_log_rho)
unified_log_rho = unified_log_rho[sort_idx]
unified_log_P = unified_log_P[sort_idx]

# Export the unified EOS to CSV file
eos_data = pd.DataFrame({
    'log_rho': unified_log_rho,
    'log_P': unified_log_P
})

# Suggested filename
filename = "data/unified_eos_magnetic_BPS-BBP-NR-NS-Polytrope_B_001.csv"
eos_data.to_csv(filename, index=False)
print(f"Exported unified EOS to {filename}")

# PLOTTING
plt.figure(figsize=(12, 8))
# plt.scatter(log_rho_bps, log_P_bps, label="BPS", c='red', s=10, alpha=0.5)
# plt.scatter(log_rho_bbp, log_P_bbp, label="BBP", c='blue', s=10, alpha=0.5)
# plt.scatter(transition_x, transition_y, label="BPS-BBP Transition", c='green', s=10)
# plt.scatter(poly_transition_x, poly_transition_y, label="BBP-Polytrope Transition", c='orange', s=10)
# plt.scatter(log_neutron_polytrope_rho, log_neutron_polytrope_P, 
#             label="Ultra Relativistic Neutron Polytrope", c='purple', s=20, marker='x'
# )
# Highlight transition regions
plt.axvspan(x_1, x_2, alpha=0.2, color='yellow', label="BPS-BBP Transition Region")
plt.axvspan(poly_x_1, poly_x_2, alpha=0.2, color='cyan', label="BBP-Polytrope Transition Region")

plt.scatter(unified_log_rho, unified_log_P, label="Final Unified EOS", c='black', s=1)

plt.xlabel(r'$\log_{10}(\rho)$ [g/cm$^3$]')
plt.ylabel(r'$\log_{10}(P)$ [dyne/cm$^2$]')
plt.title("Complete Unified EOS with Cubic Hermite Interpolation")
plt.legend()
plt.grid(True)

# plt.xlim(12, 17)
# plt.ylim(30, 38)

plt.savefig("unified_eos_magnetic_BPS-BBP-NR-NS-Polytrope_B_001.png")

# 11. Calculate and verify physical properties
# Adiabatic index
unified_drho = np.diff(unified_log_rho)
unified_dP = np.diff(unified_log_P)
unified_adiabatic_index = unified_dP / unified_drho

non_monotonic_indices = np.where(unified_adiabatic_index <= 0)[0]
if len(non_monotonic_indices) > 0:
    for idx in non_monotonic_indices:
        print(f"Non-monotonic at ρ={unified_log_rho[idx]:.4f}, Γ={unified_adiabatic_index[idx]:.6f}")

# Check monotonicity and thermodynamic stability
monotonic = np.all(unified_adiabatic_index > 0)
stable = np.all(unified_adiabatic_index >= 4/3)
print(f"Unified EOS is monotonic: {monotonic}")
print(f"Unified EOS is thermodynamically stable: {stable}")

# Plot adiabatic index
plt.figure(figsize=(10, 5))
plt.plot(unified_log_rho[:-1], unified_adiabatic_index)
plt.axvline(x_1, linestyle='--', color='red', label="BPS-BBP Transition")
plt.axvline(x_2, linestyle='--', color='red')
plt.axvline(poly_x_1, linestyle='--', color='blue', label="BBP-Polytrope Transition")
plt.axvline(poly_x_2, linestyle='--', color='blue')
plt.axhline(4/3, linestyle=':', color='green', label="Γ = 4/3 Stability Limit")
plt.xlabel(r'$\log_{10}(\rho)$ [g/cm$^3$]')
plt.ylabel(r'Adiabatic Index ($\Gamma$)')
plt.ylim(-0.5, 10.0)
plt.title("Adiabatic Index in Unified EOS")
plt.legend()
plt.grid(True)
plt.savefig("unified_eos_magnetic_BPS-BBP-NR-NS-Polytrope_B_001_adiabatic_index.png")