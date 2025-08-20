import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import brentq
from scipy.interpolate import PchipInterpolator
from scipy.integrate import cumulative_trapezoid
import matplotlib.pyplot as plt
import math

c = 2.99792458e10  # cm/s

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
bbp_data = pd.read_csv("/home/karan-kinariwala/Dropbox/KARAN/2-Areas/Education/PhD/3-Research/pulsarmhd/data/bbp.csv")
magnetic_bps_data = pd.read_csv("/home/karan-kinariwala/Dropbox/KARAN/2-Areas/Education/PhD/3-Research/pulsarmhd/data/magnetic_bps_eos_B_001.csv")
bps_data = pd.read_csv("/home/karan-kinariwala/Dropbox/KARAN/2-Areas/Education/PhD/3-Research/pulsarmhd/data/bps.csv")

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
print("CONNECTING MAGNETIC BPS AND BBP EOS's")


# ---- Sort & de-dup both tables (BPS already unique'd; keep it explicit) ----
bps_df = (
    pd.DataFrame({'x': log_rho_bps, 'y': log_P_bps})
    .dropna()
    .drop_duplicates(subset='x', keep='first')
    .sort_values('x')
)
bbp_df = (
    pd.DataFrame({'x': log_rho_bbp, 'y': log_P_bbp})
    .dropna()
    .drop_duplicates(subset='x', keep='first')
    .sort_values('x')
)

# Overlap for diagnostics (replaces the previous min/max lines)
log_rho_min = max(bps_df['x'].iloc[0],      bbp_df['x'].iloc[0])
log_rho_max = min(bps_df['x'].iloc[-1],     bbp_df['x'].iloc[-1])
print(f"BPS density range: {bps_df['x'].iloc[0]} to {bps_df['x'].iloc[-1]}")
print(f"BBP density range: {bbp_df['x'].iloc[0]} to {bbp_df['x'].iloc[-1]}")

# ---- Shape-preserving interpolants in log–log space ----
# (extrapolate=False so we don't silently use out-of-range values)
f_bps = PchipInterpolator(bps_df['x'].values, bps_df['y'].values, extrapolate=False)
f_bbp = PchipInterpolator(bbp_df['x'].values, bbp_df['y'].values, extrapolate=False)

# ---- Root-find join (same logic as Step 1, now with PCHIP) ----
a = float(log_rho_min); b = float(log_rho_max)
if not (np.isfinite(a) and np.isfinite(b) and a < b):
    raise ValueError("No valid overlap between BPS and BBP tables to perform the join.")

def g(x):  # Δ(log10 P)
    return float(f_bps(x) - f_bbp(x))

ga, gb = g(a), g(b)
# Diagnostics to understand sign pattern
xs_probe = np.linspace(a, b, 9)
g_probe = np.array([g(xx) for xx in xs_probe])
print("g(a), g(b):", ga, gb)
print("min(g) on overlap, max(g) on overlap:", np.nanmin(g_probe), np.nanmax(g_probe))

if np.isfinite(ga) and np.isfinite(gb) and ga * gb < 0:
    log_rho_join = brentq(g, a, b)      # pressure equality
    method = "root"
else:
    xs = np.linspace(a, b, 2000)
    diffs = np.abs([g(xx) for xx in xs])
    i = int(np.nanargmin(diffs))
    log_rho_join = float(xs[i])
    method = "min_distance_fallback"

print(f"Join method: {method}")
print(r'$\log_{10}(\rho)$ at join: ', log_rho_join)

log_P_bps_transition = float(f_bps(log_rho_join))
log_P_bbp_transition = float(f_bbp(log_rho_join))
print(f"log(P) from BPS at join: {log_P_bps_transition}")
print(f"log(P) from BBP at join: {log_P_bbp_transition}")
# absolute and relative mismatch for clarity
P_bps = 10.0**log_P_bps_transition
P_bbp = 10.0**log_P_bbp_transition
print(f"Pressure mismatch (linear): {P_bps - P_bbp:.6e}  (ratio BPS/BBP = {P_bps / P_bbp:.6f})")

def phi0(t):
    return 1 - (3 * (t ** 2)) + (2 * (t ** 3))

def phi1(t):
    return (3 * (t ** 2)) - (2 * (t ** 3))

def psi0(t):
    return (t ** 3) - (2 * (t ** 2)) + t

def psi1(t):
    return (t ** 3) - (t ** 2)

def cubic_hermite(x, x1, x2, f1, f2, m1, m2):
    h = x2 - x1
    t = (x - x1) / h
    return f1*phi0(t) + f2*phi1(t) + h*m1*psi0(t) + h*m2*psi1(t)

# Keep x1 near the min-|Δ| point but away from the hard boundary by a tiny margin
x0 = float(log_rho_join)
a_eps = a + 1e-4
b_eps = b - 1e-4
x_1 = max(a_eps, min(x0, b_eps))

y1 = float(f_bps(x_1))  # log10 P at start (BPS)
P1 = 10.0**y1

# target: BBP reaches slightly above P1 (e.g., +0.1%)
frac = 1e-3
P_target = P1 * (1.0 + frac)
y_target = math.log10(P_target)

# Find x2 >= x1 such that f_bbp(x2) == y_target, within BBP domain (extrapolate=False)
x2_lo = x_1
x2_hi = float(bbp_df['x'].iloc[-1])  # high end of BBP table
def h(x): 
    val = f_bbp(x)
    return float(val - y_target) if np.isfinite(val) else -1.0  # keep sign

h_lo, h_hi = h(x2_lo), h(x2_hi)

if np.isnan(h_lo):
    raise RuntimeError("BBP undefined at x1; check domain/NaNs.")

if h_hi <= 0:
    # Even at the top of BBP table we haven't reached y_target.
    # Fall back to exact equality with y1 (plateau end), which MUST exist eventually.
    # Try root for y1 instead of y_target:
    y_plateau = y1
    def h_eq(x): 
        val = f_bbp(x)
        return float(val - y_plateau) if np.isfinite(val) else -1.0
    if h_eq(x2_hi) <= 0:
        # Extreme fallback: clamp x2 to table end and we will create a flat segment (PCHIP later pieces will lift it)
        x_2 = x2_hi
        f2 = float(f_bbp(x_2))
    else:
        x_2 = brentq(h_eq, x2_lo, x2_hi)
        f2 = float(f_bbp(x_2))
else:
    x_2 = brentq(h, x2_lo, x2_hi)   # BBP reaches y_target
    f2 = float(f_bbp(x_2))

# Endpoint values
f1 = y1  # log10 P at x1 (BPS)
print(f"Blend interval (monotone target): [{x_1:.6f}, {x_2:.6f}]  (width = {x_2 - x_1:.5f} dex)")
print(f"Endpoint values: f1(BPS)={f1:.6f}, f2(BBP)={f2:.6f} (should be >= f1)")

# Endpoint slopes from PCHIP
df_bps = f_bps.derivative(); df_bbp = f_bbp.derivative()
m1_raw = float(df_bps(x_1))
m2_raw = float(df_bbp(x_2))

# --- Single-interval monotone Hermite slope limiting ---
# For increasing data (f2 >= f1), let Δ = (f2 - f1)/(x2 - x1).
# Require m1,m2 >= 0 and not too large relative to Δ to avoid overshoot.
def limit_slopes_monotone(f1, f2, x1, x2, m1, m2):
    h = x2 - x1
    if h <= 0:
        return 0.0, 0.0
    delta = (f2 - f1) / h
    if delta <= 0:  # plateau or decreasing -> choose flat to remain monotone nondecreasing
        return 0.0, 0.0
    # Clip to nonnegative
    m1 = max(0.0, m1)
    m2 = max(0.0, m2)
    # Express as α=m1/Δ, β=m2/Δ
    alpha = m1 / delta; beta = m2 / delta
    # Sufficient conditions for no overshoot on a single interval:
    # α <= 3, β <= 3 and α^2 + β^2 <= 9  (Fritsch–Carlson/Hyman-style limiter)
    # First clamp to 3
    alpha = min(alpha, 3.0); beta = min(beta, 3.0)
    s2 = alpha*alpha + beta*beta
    if s2 > 9.0:
        tau = 3.0 / math.sqrt(s2)
        alpha *= tau; beta *= tau
    return alpha*delta, beta*delta

m1, m2 = limit_slopes_monotone(f1, f2, x_1, x_2, m1_raw, m2_raw)
print(f"Endpoint slopes (raw -> limited): m1 {m1_raw:.6f} -> {m1:.6f},  m2 {m2_raw:.6f} -> {m2:.6f}")

# 3c) Build the hybrid EOS: BPS for x<x1, Hermite on [x1,x2], BBP for x>x2
# Use the original tabulated (sorted, de-duped) points to keep native resolution.
# Build hybrid arrays
bps_left_mask  = bps_df['x'].values < x_1
bbp_right_mask = bbp_df['x'].values > x_2
bps_left_x = bps_df['x'].values[bps_left_mask]
bps_left_y = bps_df['y'].values[bps_left_mask]
bbp_right_x = bbp_df['x'].values[bbp_right_mask]
bbp_right_y = bbp_df['y'].values[bbp_right_mask]

trans_x = np.linspace(x_1, x_2, 200)
trans_y = np.array([cubic_hermite(xx, x_1, x_2, f1, f2, m1, m2) for xx in trans_x])

# --- Monotone flank guard (physically dP/dρ ≥ 0) ---

# Enforce non-decreasing y on the left flank (BPS)
if bps_left_y.size:
    bps_left_y = np.maximum.accumulate(bps_left_y)  # isotonic correction
    # Ensure no downward jump into the transition start
    bps_left_y[-1] = min(bps_left_y[-1], f1)

# Enforce non-decreasing y on the right flank (BBP)
if bbp_right_y.size:
    # Ensure no downward jump out of the transition end
    bbp_right_y[0] = max(bbp_right_y[0], f2)
    bbp_right_y = np.maximum.accumulate(bbp_right_y)

hybrid_log_rho = np.concatenate([bps_left_x, trans_x, bbp_right_x])
hybrid_log_P   = np.concatenate([bps_left_y, trans_y, bbp_right_y])

# Diagnostics: monotonicity in log–log (allow tiny plateaus)
dlogP = np.diff(hybrid_log_P); dlogR = np.diff(hybrid_log_rho)
Gamma1 = dlogP / dlogR
tol = 1e-10  # tolerance in log–log slope
print(f"Gamma1 min/max on hybrid: {np.nanmin(Gamma1):.6f} / {np.nanmax(Gamma1):.6f}")
neg = np.where(Gamma1 < -tol)[0]  # strictly negative beyond tolerance is a failure
print("Monotone (non-decreasing) in log–log:", neg.size == 0)
if neg.size:
    print("First non-monotone indices:", neg[:10])

print("\nCONNECTING BBP-HYBRID AND NR NEUTRON POLYTROPE (γ=5/3)")

# --- Build PCHIP over the hybrid you just constructed ---
f_hybrid = PchipInterpolator(hybrid_log_rho, hybrid_log_P, extrapolate=False)
df_hybrid = f_hybrid.derivative()

# --- Define the NR neutron polytrope in log–log space ---
K_poly = 5.3802e9  # cgs
gamma_poly = 5.0/3.0

def logP_poly(x):           # x = log10 rho
    return math.log10(K_poly) + gamma_poly * x

def dlogP_dlogrho_poly(x):  # derivative in log–log is constant = gamma
    return gamma_poly

# Polytrope domain (in log10 rho) — keep it generous but finite
x_poly_min = 12.0
x_poly_max = 15.5

# --- Overlap between hybrid and polytrope domains ---
x_hyb_min = float(hybrid_log_rho[0])
x_hyb_max = float(hybrid_log_rho[-1])

a2 = max(x_hyb_min, x_poly_min)
b2 = min(x_hyb_max, x_poly_max)

if not (np.isfinite(a2) and np.isfinite(b2) and a2 < b2):
    raise ValueError("No valid overlap between HYBRID and POLYTROPE to perform the join.")

# --- Try pressure equality first ---
def g2(x):  # Δ(log10 P) = hybrid - polytrope
    fh = f_hybrid(x)
    return float(fh - logP_poly(x)) if np.isfinite(fh) else np.nan

g2a, g2b = g2(a2), g2(b2)

xs_probe2 = np.linspace(a2, b2, 9)
g2_probe = np.array([g2(xx) for xx in xs_probe2])
print("g2(a2), g2(b2):", g2a, g2b)
print("min(g2), max(g2) on overlap:", np.nanmin(g2_probe), np.nanmax(g2_probe))

if np.isfinite(g2a) and np.isfinite(g2b) and g2a * g2b < 0:
    x_join2 = brentq(lambda x: g2(x), a2, b2)
    method2 = "root"
else:
    xs2 = np.linspace(a2, b2, 2000)
    diffs2 = np.abs([g2(xx) for xx in xs2])
    j2 = int(np.nanargmin(diffs2))
    x_join2 = float(xs2[j2])
    method2 = "min_distance_fallback"
print(f"Join method (hybrid→poly): {method2}")
print(r'$\log_{10}(\rho)$ at join: ', x_join2)

# --- Build a short monotone blend from hybrid (left) to polytrope (right) ---
# x1 near the join on hybrid; x2 where polytrope reaches slightly above f1
x1 = max(a2 + 1e-4, min(x_join2, b2 - 1e-4))
f1J = float(f_hybrid(x1))
P1J = 10.0**f1J

frac = 1e-3             # 0.1% upward target
P2_target = P1J * (1.0 + frac)
y2_target = math.log10(P2_target)

# Find x2 ≥ x1 with logP_poly(x2) = y2_target (analytic)
x2 = (y2_target - math.log10(K_poly))/gamma_poly
if not np.isfinite(x2) or x2 <= x1:
    # Fallback: exact equality with f1J (plateau end)
    x2 = (f1J - math.log10(K_poly))/gamma_poly
    if not np.isfinite(x2) or x2 <= x1:
        # Final fallback: very short band
        x2 = min(b2, x1 + 0.03)
y2 = logP_poly(x2)

print(f"Blend interval (hybrid→poly): [{x1:.6f}, {x2:.6f}]  (width = {x2 - x1:.5f} dex)")
print(f"Endpoint values: f1(hybrid)={f1J:.6f}, f2(poly)={y2:.6f} (should be >= f1)")

# Endpoint slopes (raw), then limit for single-interval monotonicity
m1_raw = float(df_hybrid(x1))
m2_raw = dlogP_dlogrho_poly(x2)

def limit_slopes_monotone_single(f1, f2, x1, x2, m1, m2):
    h = x2 - x1
    if h <= 0:
        return 0.0, 0.0
    delta = (f2 - f1) / h
    if delta <= 0:
        return 0.0, 0.0
    m1 = max(0.0, m1)
    m2 = max(0.0, m2)
    alpha = m1 / delta; beta = m2 / delta
    alpha = min(alpha, 3.0); beta = min(beta, 3.0)
    s2 = alpha*alpha + beta*beta
    if s2 > 9.0:
        tau = 3.0 / math.sqrt(s2)
        alpha *= tau; beta *= tau
    return alpha*delta, beta*delta

m1J, m2J = limit_slopes_monotone_single(f1J, y2, x1, x2, m1_raw, m2_raw)
print(f"Endpoint slopes (raw -> limited): m1 {m1_raw:.6f} -> {m1J:.6f},  m2 {m2_raw:.6f} -> {m2J:.6f}")

# Hermite over [x1,x2]
trans2_x = np.linspace(x1, x2, 200)
trans2_y = np.array([cubic_hermite(xx, x1, x2, f1J, y2, m1J, m2J) for xx in trans2_x])

# Left of x1: take hybrid; Right of x2: take polytrope
hyb_left_mask  = hybrid_log_rho < x1
poly_right_mask = np.array([True,])  # we will generate poly points directly

# Generate a polytrope tail (log grid) starting just beyond x2 up to x_poly_max
poly_tail_x = np.linspace(x2, x_poly_max, 400)
poly_tail_y = logP_poly(poly_tail_x)

# Monotone flank guards
hyb_left_x = hybrid_log_rho[hyb_left_mask]
hyb_left_y = hybrid_log_P[hyb_left_mask]
if hyb_left_y.size:
    hyb_left_y = np.maximum.accumulate(hyb_left_y)
    hyb_left_y[-1] = min(hyb_left_y[-1], f1J)

poly_tail_y[0] = max(poly_tail_y[0], trans2_y[-1])
poly_tail_y = np.maximum.accumulate(poly_tail_y)

# Final unified EOS (log space)
unified_log_rho = np.concatenate([hyb_left_x, trans2_x, poly_tail_x[1:]])  # avoid duplicate at x2
unified_log_P   = np.concatenate([hyb_left_y, trans2_y, poly_tail_y[1:]])

# Diagnostics
dlogP = np.diff(unified_log_P); dlogR = np.diff(unified_log_rho)
Gamma1 = dlogP / dlogR
tol = 1e-10
print(f"Gamma1 min/max (final): {np.nanmin(Gamma1):.6f} / {np.nanmax(Gamma1):.6f}")
neg = np.where(Gamma1 < -tol)[0]
print("Final EOS monotone (non-decreasing) in log–log:", neg.size == 0)
if neg.size:
    print("First non-monotone indices (final):", neg[:10])


# ---- STEP 5: Export solver-ready EOS ----
# Convert back to linear units for TOV and compute Gamma1 and dP/drho

rho = 10.0**unified_log_rho            # g cm^-3
P   = 10.0**unified_log_P              # dyne cm^-2

# Gamma1 = d ln P / d ln rho (base independent)
# Use central differences where possible; fall back to one-sided at ends
Gamma1 = np.empty_like(P)
Gamma1[1:-1] = (unified_log_P[2:] - unified_log_P[:-2]) / (unified_log_rho[2:] - unified_log_rho[:-2])
Gamma1[0]    = (unified_log_P[1] - unified_log_P[0]) / (unified_log_rho[1] - unified_log_rho[0])
Gamma1[-1]   = (unified_log_P[-1] - unified_log_P[-2]) / (unified_log_rho[-1] - unified_log_rho[-2])

# dP/drho = Gamma1 * P / rho  (derivation: dlnP = Gamma1 dlnrho)
dP_drho = Gamma1 * (P / rho)

# Optional: enforce strictly increasing P(ρ) for codes that demand it
strict = False  # set True if your downstream needs strictly increasing P
if strict:
    # tiny seam nudges to avoid any zero finite-difference
    eps_logP = 1e-7
    # identify the two seam x's used above: x_1 (BPS->BBP) and x2 (hybrid->poly)
    # nudge immediate neighbors if needed
    for seam_x in [x_1, x2]:  # x2 from the hybrid->poly block
        i = np.searchsorted(unified_log_rho, seam_x)
        i = min(max(i, 1), len(unified_log_rho)-1)
        if unified_log_P[i] <= unified_log_P[i-1]:
            unified_log_P[i] = unified_log_P[i-1] + eps_logP
    # recompute linear arrays and derivatives after nudges
    P   = 10.0**unified_log_P
    dP_drho = Gamma1 * (P / rho)

# Energy density ε(ρ) via the cold-EOS identity
I = cumulative_trapezoid(P / rho**2, rho, initial=0.0)             # ∫ P/ρ^2 dρ
epsilon = rho * (c**2 + I)                                 # erg cm^-3

# Sound speed and causality check: c_s^2 = dP/dε
cs2 = dP_drho / ((epsilon + P) / rho)                      # dimensionless (fraction of c^2)
print(f"Causality OK: {np.all(cs2 <= 1.0)}  (max c_s^2 = {cs2.max():.4f})")

# Optional: region flags for debugging/plots
region = np.full_like(unified_log_rho, fill_value="poly", dtype=object)
region[unified_log_rho < x1] = "hybrid"
region[(unified_log_rho >= x1) & (unified_log_rho <= x2)] = "hyb_poly_transition"
region[(unified_log_rho < x_2)] = "bps_bbp_or_transition"  # coarse label for left half

# Export a TOV-ready table
out = pd.DataFrame({
    "rho": rho,               # g/cm^3
    "P": P,                   # dyne/cm^2
    "epsilon": epsilon,       # erg/cm^3
    "Gamma1": Gamma1,
    "dP_drho": dP_drho,
    "cs2": cs2                # in units of c^2
})

# Basic integrity checks
assert np.all(np.diff(P) >= -1e-20), "Non-monotone P detected."
assert np.all(Gamma1 >= -1e-12), "Negative Gamma1 detected."

# Write
out_path = "/home/karan-kinariwala/Dropbox/KARAN/2-Areas/Education/PhD/3-Research/pulsarmhd/data/unified_eos_magnetic_BPS-BBP-Polytrope_B_001.csv"
out.to_csv(out_path, index=False)
print(f"Exported unified EOS to {out_path}  (rows={len(out)})")


# ===================== PLOTTING =====================

# --- 1) EOS curve: log10 P vs log10 rho with blend bands ---
plt.figure(figsize=(10, 6))

# Raw source tables (lines to show overall shape)
plt.plot(bps_df['x'].values, bps_df['y'].values, lw=1.0, alpha=0.6, label="BPS (magnetic)")
plt.plot(bbp_df['x'].values, bbp_df['y'].values, lw=1.0, alpha=0.6, label="BBP")

# Hybrid segment (after BPS→BBP seam) and final unified EOS
plt.plot(hybrid_log_rho, hybrid_log_P, lw=1.2, alpha=0.9, label="Hybrid (BPS→BBP)")
plt.plot(unified_log_rho, unified_log_P, lw=2.0, label="Unified EOS")

# Highlight transition regions
plt.axvspan(x_1, x_2, alpha=0.15, color='yellow', label="BPS→BBP blend")
plt.axvspan(x1,  x2,  alpha=0.15, color='cyan',   label="Hybrid→Poly blend")

plt.xlabel(r'$\log_{10}(\rho)\ \,[\mathrm{g\,cm^{-3}}]$')
plt.ylabel(r'$\log_{10}(P)\ \,[\mathrm{dyn\,cm^{-2}}]$')
plt.title("Unified EOS (log–log)")
plt.grid(True, alpha=0.3)
plt.legend(loc="best", frameon=True)
plt.tight_layout()
plt.savefig("/home/karan-kinariwala/Dropbox/KARAN/2-Areas/Education/PhD/3-Research/pulsarmhd/plots/unified_eos_magnetic_BPS-BBP-Polytrope_B_001.png", dpi=200)


# --- 2) Adiabatic index Γ1 = d ln P / d ln rho ---
# Use a stable derivative; if you've already computed `Gamma1`, reuse it.
Gamma1_plot = np.gradient(unified_log_P, unified_log_rho)

plt.figure(figsize=(10, 5))
plt.plot(unified_log_rho, Gamma1_plot, lw=1.6, label=r'$\Gamma_1$')

# Mark blend regions
plt.axvspan(x_1, x_2, alpha=0.15, color='yellow', label="BPS→BBP blend")
plt.axvspan(x1,  x2,  alpha=0.15, color='cyan',   label="Hybrid→Poly blend")

# Reference: Γ = 4/3 (orientation only; not a local stability criterion)
plt.axhline(4.0/3.0, linestyle=':', label=r'$\Gamma=4/3$ (reference)')

# Y-limits tuned to suppress rare large spikes but retain structure
ymin = max(0.0, np.nanpercentile(Gamma1_plot, 0.5) - 0.2)
ymax = np.nanpercentile(Gamma1_plot, 99.5)
ymax = float(np.clip(ymax, 2.0, 6.0))  # cap for readability; adjust if needed
plt.ylim(ymin, ymax)

plt.xlabel(r'$\log_{10}(\rho)\ \,[\mathrm{g\,cm^{-3}}]$')
plt.ylabel(r'$\Gamma_1 = d\ln P / d\ln \rho$')
plt.title("Adiabatic Index of Unified EOS")
plt.grid(True, alpha=0.3)
plt.legend(loc="best", frameon=True)
plt.tight_layout()
plt.savefig("/home/karan-kinariwala/Dropbox/KARAN/2-Areas/Education/PhD/3-Research/pulsarmhd/plots/unified_eos_magnetic_BPS-BBP-Polytrope_B_001_gamma.png", dpi=200)