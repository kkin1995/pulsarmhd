import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator
import matplotlib.pyplot as plt

speed_of_light = 29979245800 # cm s-1

def phi0(t): return 1 - 3*t**2 + 2*t**3
def phi1(t): return 3*t**2 - 2*t**3
def psi0(t): return t**3 - 2*t**2 + t
def psi1(t): return t**3 - t**2

def cubic_hermite(x, x1, x2, f1, f2, m1, m2):
    """Scalar cubic Hermite in log–log space on [x1,x2]."""
    h = x2 - x1
    t = (x - x1) / h
    return f1*phi0(t) + f2*phi1(t) + h*m1*psi0(t) + h*m2*psi1(t)

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

sahu_basu_datta_eos = pd.read_csv(
    "/home/karan-kinariwala/Dropbox/KARAN/2-Areas/Education/PhD/3-Research/pulsarmhd/data/sahu-basu-datta.dat",
    sep=r'\s+', # whitespace character
    header=None,
    dtype=float
)

print(f"Number of NaN Values in Sahu Basu Datta: {sahu_basu_datta_eos.isna().sum(axis = 0)}")

magnetic_bps_eos = pd.read_csv(
    "/home/karan-kinariwala/Dropbox/KARAN/2-Areas/Education/PhD/3-Research/pulsarmhd/data/bps_single_b_b_0e_00.csv",
    header=0
)

print(f"Number of NaN Values in Magnetic BPS: {magnetic_bps_eos.isna().sum(axis = 0)}")

bbp_eos = pd.read_csv(
    "/home/karan-kinariwala/Dropbox/KARAN/2-Areas/Education/PhD/3-Research/pulsarmhd/data/bbp.csv",
    header = 0
)

print(f"Number of NaN Values in BBP: {bbp_eos.isna().sum(axis = 0)}")

neutron_drip = 4.46e11
sbd_mass_density_min = sahu_basu_datta_eos.iloc[:, 1].min() * 1e14
bbp_mass_density_max = bbp_eos.iloc[:, 0].max()

magnetic_bps_eos = magnetic_bps_eos[magnetic_bps_eos.iloc[:, 1] <= np.log10(neutron_drip)]
bbp_eos = bbp_eos[(bbp_eos.iloc[:, 0] >= neutron_drip) & (bbp_eos.iloc[:, 0] <= bbp_mass_density_max)]
sahu_basu_datta_eos = sahu_basu_datta_eos[sahu_basu_datta_eos.iloc[:, 1] >= (bbp_mass_density_max / 1e14)]

sbd_number_density = sahu_basu_datta_eos.iloc[:, 0].values * 1e14
sbd_mass_density = sahu_basu_datta_eos.iloc[:, 1].values * 1e14
sbd_pressure = sahu_basu_datta_eos.iloc[:, 2].values * (speed_of_light ** 2) * 1e14

print(f"Density Range (Sahu Basu Datta): Min: {min(sbd_mass_density)} g cm-3 | Max: {max(sbd_mass_density)} g cm-3")

sbd_log_n = np.log10(sbd_number_density)
sbd_log_rho = np.log10(sbd_mass_density)
sbd_log_P = np.log10(sbd_pressure)

magnetic_bps_log_n = magnetic_bps_eos.iloc[:, 0].values
magnetic_bps_log_rho = magnetic_bps_eos.iloc[:, 1].values
magnetic_bps_log_P = magnetic_bps_eos.iloc[:, 2].values

print(f"Density Range (Magnetic BPS): Min: {min(10 ** magnetic_bps_log_rho)} g cm-3 | Max: {max(10 ** magnetic_bps_log_rho)} g cm-3")

# rho,P,nb,A,Z,Gamma
bbp_mass_density = bbp_eos.iloc[:, 0].values
bbp_pressure = bbp_eos.iloc[:, 1].values

bbp_log_rho = np.log10(bbp_mass_density)
bbp_log_P = np.log10(bbp_pressure)

print(f"Density Range (BBP): Min: {min(bbp_mass_density)} g cm-3 | Max: {max(bbp_mass_density)} g cm-3")

f_magnetic_bps = PchipInterpolator(magnetic_bps_log_rho, magnetic_bps_log_P, extrapolate=False)
f_bbp = PchipInterpolator(bbp_log_rho, bbp_log_P, extrapolate=False)

composite_log_rho = np.concatenate([sbd_log_rho, bbp_log_rho, magnetic_bps_log_rho])
composite_log_P = np.concatenate([sbd_log_P, bbp_log_P, magnetic_bps_log_P])

# At Magnetic BPS → BBP transition
mbps_max_idx = np.argmax(magnetic_bps_log_rho)
bbp_min_idx = np.argmin(bbp_log_rho)
print(f"Magnetic BPS max: log(ρ)={magnetic_bps_log_rho[mbps_max_idx]:.2f}, log(P)={magnetic_bps_log_P[mbps_max_idx]:.2f}")
print(f"BBP min: log(ρ)={bbp_log_rho[bbp_min_idx]:.2f}, log(P)={bbp_log_P[bbp_min_idx]:.2f}")

# At BBP → SBD transition
bbp_max_idx = np.argmax(bbp_log_rho)
sbd_min_idx = np.argmin(sbd_log_rho)
print(f"BBP max: log(ρ)={bbp_log_rho[bbp_max_idx]:.2f}, log(P)={bbp_log_P[bbp_max_idx]:.2f}")
print(f"SBD min: log(ρ)={sbd_log_rho[sbd_min_idx]:.2f}, log(P)={sbd_log_P[sbd_min_idx]:.2f}")

composite_eos_df = pd.DataFrame(
    {
        "log_rho": composite_log_rho,
        "log_P": composite_log_P
    }
)
composite_eos_df.to_csv("data/sahu_basu_datta_bbp_magnetic_bps_b_0e_00.csv", index=False)

fig, ax = plt.subplots(figsize=(10, 5))
# ax.scatter(composite_log_rho, composite_log_P, s = 5, c = 'blue')
ax.scatter(sbd_log_rho, sbd_log_P, s = 5, c = 'red', label = "Sahu - Basu - Datta")
ax.scatter(bbp_log_rho, bbp_log_P, s = 5, c = 'green', label = "Baym - Bethe - Pethick")
ax.scatter(magnetic_bps_log_rho, magnetic_bps_log_P, s = 5, c = 'purple', label = "Magnetic Baym - Pethick - Sutherland")
ax.set_xlabel(r'$\log \rho$')
ax.set_ylabel(r'$\log P$')
ax.set_title("Composite EOS (SBD-BBP-Magnetic-BPS)")
plt.legend()
plt.show()
