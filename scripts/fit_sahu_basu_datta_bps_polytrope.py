import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator
from scipy.integrate import cumulative_trapezoid
import matplotlib.pyplot as plt
import math

speed_of_light = 29979245800  # cm s-1
c2 = speed_of_light**2  # (cm/s)^2

def clean_monotone_xy(x, y):
    """Return finite, strictly increasing x and aligned y (first occurrence kept on duplicates)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    finite = np.isfinite(x) & np.isfinite(y)
    x = x[finite]; y = y[finite]
    if x.size == 0:
        raise ValueError("No finite points remain after cleaning.")
    order = np.argsort(x)
    x = x[order]; y = y[order]
    if x.size >= 2:
        keep = np.ones_like(x, dtype=bool)
        keep[1:] = x[1:] > x[:-1]  # strictly increasing
        x = x[keep]; y = y[keep]
    if x.size < 2:
        raise ValueError("Need at least two strictly increasing x values after cleaning.")
    return x, y

def phi0(t): return 1 - 3*t**2 + 2*t**3
def phi1(t): return 3*t**2 - 2*t**3
def psi0(t): return t**3 - 2*t**2 + t
def psi1(t): return t**3 - t**2

def cubic_hermite(x, x1, x2, f1, f2, m1, m2):
    """Scalar cubic Hermite in log–log space on [x1,x2]. m's are dy/dx."""
    h = x2 - x1
    t = (x - x1) / h
    return f1*phi0(t) + f2*phi1(t) + h*m1*psi0(t) + h*m2*psi1(t)

def fritsch_carlson_monotonic_slopes(x1, x2, y1, y2, m1, m2):
    dx = x2 - x1
    dy = y2 - y1
    if abs(dx) < 1e-10: return 0.0, 0.0
    delta = dy / dx
    if abs(delta) < 1e-10: return 0.0, 0.0
    alpha = m1 / delta; beta = m2 / delta
    if alpha < 0 or beta < 0: return 0.0, 0.0
    sum_ab = alpha + beta
    cond2 = alpha + 2*beta - 3
    cond3 = 2*alpha + beta - 3
    if sum_ab == 2: return m1, m2
    cond1 = False
    if sum_ab > 2:
        term1 = 2*alpha + beta - 3
        term2 = term1*term1/(3*(sum_ab - 2))
        cond1 = alpha*term2 > 0
    if cond1 or cond2 <= 0 or cond3 <= 0: return m1, m2
    if sum_ab > 3:
        tau = 3.0 / sum_ab
        return tau*m1, tau*m2
    ab2 = alpha*alpha + beta*beta
    if ab2 > 9:
        tau = 3.0 / math.sqrt(ab2)
        return tau*m1, tau*m2
    return 0.9*m1, 0.9*m2

def build_c1_bridge(xL, yL, dL, xR, yR, dR, n_mid=2):
    """
    Create n_mid interior points on a C¹ cubic Hermite between (xL,yL) and (xR,yR).
    Slopes dL,dR are dy/dx; clip with Fritsch–Carlson to preserve monotonicity.
    """
    if not np.isfinite([xL,yL,dL,xR,yR,dR]).all():
        raise ValueError("Non-finite endpoint or slope in bridge.")
    if not (xR > xL and yR >= yL):
        raise ValueError("Bridge endpoints not monotone in x or y.")
    mL, mR = fritsch_carlson_monotonic_slopes(xL, xR, yL, yR, dL, dR)
    xs = np.linspace(xL, xR, n_mid + 2)[1:-1]  # strictly interior
    ys = np.array([cubic_hermite(x, xL, xR, yL, yR, mL, mR) for x in xs])
    return xs, ys

print("=" * 70)
print("  Composite EOS Builder: SBD + BBP + Magnetic BPS")
print("=" * 70)

# ---------- Load EOS pieces ----------
sahu_basu_datta_eos = pd.read_csv(
    "/home/karan-kinariwala/Dropbox/KARAN/2-Areas/Education/PhD/3-Research/pulsarmhd/data/sahu-basu-datta.dat",
    sep=r'\s+', header=None, dtype=float
)
# Check for data quality (silent unless issues found)
n_nan_sbd = sahu_basu_datta_eos.isna().sum().sum()
if n_nan_sbd > 0:
    print(f"[warn] Sahu-Basu-Datta has {n_nan_sbd} NaN values")

# Magnetic field parameter (used for filenames and plot labels)
B_field_str = "1e-03"  # Change this to match your magnetic BPS file
outfile = f"data/sahu_basu_datta_bbp_magnetic_bps_b_{B_field_str}.csv"

magnetic_bps_eos = pd.read_csv(
    f"/home/karan-kinariwala/Dropbox/KARAN/2-Areas/Education/PhD/3-Research/pulsarmhd/data/bps_single_b_b_{B_field_str}.csv",
    header=0
)
n_nan_bps = magnetic_bps_eos.isna().sum().sum()
if n_nan_bps > 0:
    print(f"[warn] Magnetic BPS has {n_nan_bps} NaN values")

# Attempt to obtain epsilon from BPS (optional anchor)
bps_eps_raw = magnetic_bps_eos['eps'].to_numpy(dtype=float) if 'eps' in magnetic_bps_eos.columns else None

bbp_eos = pd.read_csv(
    "/home/karan-kinariwala/Dropbox/KARAN/2-Areas/Education/PhD/3-Research/pulsarmhd/data/bbp.csv",
    header = 0
)
n_nan_bbp = bbp_eos.isna().sum().sum()
if n_nan_bbp > 0:
    print(f"[warn] BBP has {n_nan_bbp} NaN values")

# ---------- Splicing windows ----------
neutron_drip = 4.46e11
sbd_mass_density_min = sahu_basu_datta_eos.iloc[:, 1].min() * 1e14
bbp_mass_density_max = bbp_eos.iloc[:, 0].max()

# BPS: handle new format - compute logs from linear columns if needed
if 'log_rho' not in magnetic_bps_eos.columns:
    print("Detected new BPS format (linear columns), computing logs...")
    # Optional: keep only converged rows if the column exists and is usable
    if 'converged' in magnetic_bps_eos.columns and magnetic_bps_eos['converged'].notna().any():
        magnetic_bps_eos = magnetic_bps_eos[magnetic_bps_eos['converged'] == 1]

    # Require positive rho/P; include nB only if present
    need_cols = {'rho', 'P'}
    if not need_cols.issubset(magnetic_bps_eos.columns):
        raise RuntimeError(f"BPS file missing required columns {need_cols}")

    pos_mask = (magnetic_bps_eos['rho'] > 0) & (magnetic_bps_eos['P'] > 0)
    if 'nB' in magnetic_bps_eos.columns:
        pos_mask &= (magnetic_bps_eos['nB'] > 0)

    magnetic_bps_eos = magnetic_bps_eos[pos_mask].copy()
    magnetic_bps_eos['log_rho'] = np.log10(magnetic_bps_eos['rho'].astype(float))
    magnetic_bps_eos['log_P']   = np.log10(magnetic_bps_eos['P'].astype(float))
    if 'nB' in magnetic_bps_eos.columns:
        magnetic_bps_eos['log_n'] = np.log10(magnetic_bps_eos['nB'].astype(float))

    # Diagnostic: check density range before filtering (only if verbose)
    # print(f"BPS density range: {magnetic_bps_eos['rho'].min():.2e} to {magnetic_bps_eos['rho'].max():.2e} g/cm^3")

# Replace ±inf with NaN and drop
for col in ['log_rho', 'log_P']:
    if col in magnetic_bps_eos.columns:
        magnetic_bps_eos[col] = magnetic_bps_eos[col].replace([np.inf, -np.inf], np.nan)
magnetic_bps_eos = magnetic_bps_eos.dropna(subset=['log_rho','log_P'])
magnetic_bps_eos = magnetic_bps_eos[magnetic_bps_eos['log_rho'].to_numpy() <= np.log10(neutron_drip)]

# Safety check: ensure we have data after filtering
if len(magnetic_bps_eos) == 0:
    raise RuntimeError(f"No valid magnetic BPS data after filtering. Check input file and ensure it contains "
                       f"finite, positive values for rho, P, and nB up to neutron drip density {neutron_drip:.2e} g/cm^3")

# Align epsilon to filtered BPS data (if present)
if bps_eps_raw is not None and 'eps' in magnetic_bps_eos.columns:
    bps_eps_aligned = magnetic_bps_eos['eps'].to_numpy(dtype=float)
    bps_lr_src = magnetic_bps_eos['log_rho'].to_numpy(dtype=float)
else:
    bps_eps_aligned = None
    bps_lr_src = None

# BBP: keep in [drip, bbp_max]
bbp_eos = bbp_eos[(bbp_eos.iloc[:, 0] >= neutron_drip) & (bbp_eos.iloc[:, 0] <= bbp_mass_density_max)]

# SBD: keep above BBP max (its table is in units of 1e14)
sahu_basu_datta_eos = sahu_basu_datta_eos[sahu_basu_datta_eos.iloc[:, 1] >= (bbp_mass_density_max / 1e14)]

# ---------- Transform to logs ----------
sbd_number_density = sahu_basu_datta_eos.iloc[:, 0].values * 1e14
sbd_mass_density   = sahu_basu_datta_eos.iloc[:, 1].values * 1e14
sbd_pressure       = sahu_basu_datta_eos.iloc[:, 2].values * (speed_of_light ** 2) * 1e14

sbd_log_n   = np.log10(sbd_number_density)
sbd_log_rho = np.log10(sbd_mass_density)
sbd_log_P   = np.log10(sbd_pressure)

magnetic_bps_log_rho = magnetic_bps_eos['log_rho'].values
magnetic_bps_log_P   = magnetic_bps_eos['log_P'].values

# BBP to logs (guard positivity)
bbp_mass_density = bbp_eos.iloc[:, 0].values
bbp_pressure     = bbp_eos.iloc[:, 1].values
mask_bbp = np.isfinite(bbp_mass_density) & np.isfinite(bbp_pressure) & (bbp_mass_density > 0) & (bbp_pressure > 0)
bbp_mass_density = bbp_mass_density[mask_bbp]
bbp_pressure     = bbp_pressure[mask_bbp]
bbp_log_rho = np.log10(bbp_mass_density)
bbp_log_P   = np.log10(bbp_pressure)

# ---------- Minimal cleaning for PCHIP ----------
magnetic_bps_log_rho, magnetic_bps_log_P = clean_monotone_xy(magnetic_bps_log_rho, magnetic_bps_log_P)

# Capture reference epsilon from last kept BPS point (after cleaning)
eps_ref_from_bps = np.nan
if bps_eps_aligned is not None and bps_lr_src is not None and magnetic_bps_log_rho.size > 0:
    # Take epsilon at the nearest original BPS log_rho to the last kept BPS log_rho
    ref_lr = magnetic_bps_log_rho[-1]
    idx_ref = int(np.nanargmin(np.abs(bps_lr_src - ref_lr)))
    val = bps_eps_aligned[idx_ref] if 0 <= idx_ref < bps_eps_aligned.size else np.nan
    eps_ref_from_bps = val if (np.isfinite(val) and val > 0.0) else np.nan

bbp_log_rho,          bbp_log_P          = clean_monotone_xy(bbp_log_rho,          bbp_log_P)
# *** KEY FIX: clean SBD too (duplicates/out-of-order lines cause your error) ***
sbd_log_rho,          sbd_log_P          = clean_monotone_xy(sbd_log_rho,          sbd_log_P)

# ---------- Trim BBP leading rows to avoid pressure drop at BPS→BBP ----------
last_bps_rho = magnetic_bps_log_rho[-1]
last_bps_P   = magnetic_bps_log_P[-1]
tol = 1e-9
i0 = np.searchsorted(bbp_log_rho, last_bps_rho - tol, side='left')
while i0 < bbp_log_rho.size and bbp_log_P[i0] + tol < last_bps_P:
    i0 += 1
if i0 >= bbp_log_rho.size:
    raise RuntimeError("After trimming, BBP segment is empty. Check tables/splice ranges.")
if i0 > 0:
    bbp_log_rho = bbp_log_rho[i0:]
    bbp_log_P   = bbp_log_P[i0:]

# ---------- Ensure BBP→SBD x-separation and no pressure drop ----------
# If SBD starts at or before last BBP x, advance SBD start; also avoid P drop.
last_bbp_rho = bbp_log_rho[-1]
last_bbp_P   = bbp_log_P[-1]
j0 = 0
while j0 < sbd_log_rho.size and (sbd_log_rho[j0] <= last_bbp_rho + tol or sbd_log_P[j0] + tol < last_bbp_P):
    j0 += 1
if j0 >= sbd_log_rho.size:
    raise RuntimeError("After trimming, SBD segment is empty. Check tables/splice ranges.")
if j0 > 0:
    sbd_log_rho = sbd_log_rho[j0:]
    sbd_log_P   = sbd_log_P[j0:]

# ---------- Build PCHIPs and their derivatives (now safe) ----------
f_magnetic_bps = PchipInterpolator(magnetic_bps_log_rho, magnetic_bps_log_P, extrapolate=False)
f_bbp          = PchipInterpolator(bbp_log_rho,          bbp_log_P,          extrapolate=False)
f_sbd          = PchipInterpolator(sbd_log_rho,          sbd_log_P,          extrapolate=False)
df_magnetic_bps = f_magnetic_bps.derivative()
df_bbp          = f_bbp.derivative()
df_sbd          = f_sbd.derivative()

# ---------- C¹ bridges at both joins ----------
# BPS → BBP
xL, yL = magnetic_bps_log_rho[-1], magnetic_bps_log_P[-1]
xR, yR = bbp_log_rho[0],            bbp_log_P[0]
if not (yR >= yL):  # safety (should be true after BBP trim)
    raise RuntimeError("Non-monotone join BPS→BBP (P decreased after trim).")
dL = float(df_magnetic_bps(xL))
dR = float(df_bbp(xR))
xb1, yb1 = build_c1_bridge(xL, yL, dL, xR, yR, dR, n_mid=2)

# BBP → SBD
xL2, yL2 = bbp_log_rho[-1], bbp_log_P[-1]
xR2, yR2 = sbd_log_rho[0],  sbd_log_P[0]
if not (xR2 > xL2 and yR2 >= yL2):
    raise RuntimeError("Non-monotone join BBP→SBD (check SBD trim).")
dL2 = float(df_bbp(xL2))
dR2 = float(df_sbd(xR2))
xb2, yb2 = build_c1_bridge(xL2, yL2, dL2, xR2, yR2, dR2, n_mid=2)

# ---------- Composite EOS (include bridge points) ----------
comp_log_rho = np.concatenate([magnetic_bps_log_rho, xb1, bbp_log_rho, xb2, sbd_log_rho])
comp_log_P   = np.concatenate([magnetic_bps_log_P,   yb1, bbp_log_P,   yb2, sbd_log_P])

# Final tidy: strictly increasing and finite
comp_log_rho, comp_log_P = clean_monotone_xy(comp_log_rho, comp_log_P)

# Post-bridge monotonicity enforcement: use cumulative maximum to ensure strict monotonicity
# This can happen when cubic Hermite creates local maxima at bridge points
P_lin_temp = 10.0**comp_log_P
P_original = P_lin_temp.copy()

# Apply cumulative maximum: P can never decrease
P_lin_temp = np.maximum.accumulate(P_lin_temp)

# Report corrections
corrections = np.where(P_lin_temp != P_original)[0]
if len(corrections) > 0:
    for i in corrections[:3]:  # Show first 3
        print(f"[bridge-fix] Monotonicity correction at i={i}, logρ={comp_log_rho[i]:.6f}: "
              f"{P_original[i]:.4e} → {P_lin_temp[i]:.4e}")
    print(f"[bridge-fix] Applied {len(corrections)} monotonicity corrections")

comp_log_P = np.log10(P_lin_temp)

# ---------- Build epsilon consistently from composite P(rho), anchored at BPS if available ----------
rho_lin = 10.0**comp_log_rho         # g/cm^3
P_lin   = 10.0**comp_log_P           # erg/cm^3
integrand = P_lin / (rho_lin**2)     # P / rho^2

# Cumulative integral ∫ P/ρ² dρ
I_from_start = cumulative_trapezoid(integrand, rho_lin, initial=0.0)

# Find reference point (last kept BPS point in composite)
ref_idx = 0
if magnetic_bps_log_rho.size > 0:
    ref_lr = magnetic_bps_log_rho[-1]
    ref_idx = int(np.nanargmin(np.abs(comp_log_rho - ref_lr)))

rho_ref = rho_lin[ref_idx]

# Determine integration constant C
if np.isfinite(eps_ref_from_bps):
    C = (eps_ref_from_bps / rho_ref) - I_from_start[ref_idx]
else:
    C = c2 - I_from_start[0]  # Cold fallback: ε/ρ = c² at first point

epsilon = rho_lin * (C + I_from_start)  # erg/cm^3

# Gentle floor to prevent tiny negative excursions from floating-point noise (cold anchor only)
if not np.isfinite(eps_ref_from_bps):
    eps_floor = rho_lin * c2
    epsilon = np.maximum(eps_floor * (1 - 1e-12), epsilon)

# Diagnostic output
print(f"[epsilon] anchor = {'BPS eps' if np.isfinite(eps_ref_from_bps) else 'cold (rho c^2)'}")
print(f"[epsilon] ref index {ref_idx}, rho_ref={rho_ref:.4e} g/cm^3")

# ---------- Join sanity prints (quick continuity check) ----------
def nearest_idx(arr, x):
    return int(np.nanargmin(np.abs(arr - x)))

i_bpsL = nearest_idx(comp_log_rho, magnetic_bps_log_rho[-1])
i_bbpR = nearest_idx(comp_log_rho, bbp_log_rho[0])
i_sbdR = nearest_idx(comp_log_rho, sbd_log_rho[0])

def dump(i, tag):
    print(f"[{tag}] logρ={comp_log_rho[i]:.6f}  P={10**comp_log_P[i]:.6e}  "
          f"ε={epsilon[i]:.6e}  ε/ρc²={epsilon[i]/(rho_lin[i]*c2):.6e}")

print("\n=== Join Continuity Checks ===")
dump(i_bpsL, "BPS→BBP (left) ")
dump(i_bbpR, "BPS→BBP (right)")
dump(i_sbdR, "BBP→SBD (right)")

# ---------- Post-build validation (cheap physics checks) ----------
print("\n=== Physics Validation ===")

# Monotonicity in linear space
drho = np.diff(rho_lin)
dP_check = np.diff(P_lin)
deps = np.diff(epsilon)

rho_mono = np.all(drho > 0)
P_mono = np.all(dP_check >= 0)
eps_mono = np.all(deps >= 0)

if not (rho_mono and P_mono and eps_mono):
    # Find where monotonicity breaks
    if not rho_mono:
        bad_idx = np.where(drho <= 0)[0]
        print(f"[validate] ✗ Non-monotone ρ at indices: {bad_idx[:5]} (showing first 5)")
    if not P_mono:
        bad_idx = np.where(dP_check < 0)[0]
        print(f"[validate] ✗ Non-monotone P at indices: {bad_idx[:5]} (showing first 5)")
        for i in bad_idx[:3]:
            print(f"  i={i}: logρ={comp_log_rho[i]:.6f}, P[i]={P_lin[i]:.6e}, P[i+1]={P_lin[i+1]:.6e}, dP={dP_check[i]:.6e}")
    if not eps_mono:
        bad_idx = np.where(deps < 0)[0]
        print(f"[validate] ✗ Non-monotone ε at indices: {bad_idx[:5]} (showing first 5)")
        for i in bad_idx[:3]:
            print(f"  i={i}: logρ={comp_log_rho[i]:.6f}, ε[i]={epsilon[i]:.6e}, ε[i+1]={epsilon[i+1]:.6e}, dε={deps[i]:.6e}")
    raise RuntimeError("[validate] Non-monotone rho/P/epsilon after splice; check inputs.")
print("[validate] ✓ Monotonicity: ρ, P, ε all strictly/weakly increasing")

# Causality/stability (finite-diff estimates)
dP = np.diff(P_lin)
dE = np.diff(epsilon)
cs2 = dP / dE  # ≈ dP/dε (sound speed squared in natural units)

if np.any(cs2 <= 0):
    print("[warn] Found non-positive sound speed squared somewhere (dP/dε ≤ 0).")
else:
    print(f"[validate] ✓ Stability: dP/dε > 0 everywhere (min cs² = {cs2.min():.6f})")

if np.any(cs2 >= 1.0):
    print(f"[warn] Causality breach somewhere (dP/dε ≥ 1). Inspect joins/inputs. Max cs² = {cs2.max():.6f}")
else:
    print(f"[validate] ✓ Causality: dP/dε < 1 everywhere (max cs² = {cs2.max():.6f})")

pd.DataFrame({"log_rho": comp_log_rho, "log_P": comp_log_P, "epsilon": epsilon}).to_csv(outfile, index=False)

# ---------- Summary ----------
print(f"\n{'='*70}")
print(f"  Output: {outfile}")
print(f"  Total points: {len(comp_log_rho)} (BPS: {len(magnetic_bps_log_rho)}, BBP: {len(bbp_log_rho)}, SBD: {len(sbd_log_rho)})")
print(f"  Density range: 10^{comp_log_rho[0]:.2f} to 10^{comp_log_rho[-1]:.2f} g/cm³")
print(f"{'='*70}")

# ---------- Plot ----------
fig, ax = plt.subplots(figsize=(10, 5))
ax.scatter(sbd_log_rho, sbd_log_P, s=5, c='red',    label="Sahu–Basu–Datta (core, trimmed)")
ax.scatter(bbp_log_rho, bbp_log_P, s=5, c='green',  label="Baym–Bethe–Pethick (inner crust, trimmed)")
ax.scatter(magnetic_bps_log_rho, magnetic_bps_log_P, s=5, c='purple', label=f"Magnetic BPS (outer crust, b={B_field_str})")
ax.plot(xb1, yb1, 'k.', ms=6, label='C¹ bridge (BPS→BBP)')
ax.plot(xb2, yb2, 'k.', ms=6, label='C¹ bridge (BBP→SBD)')
ax.set_xlabel(r'$\log_{10}\rho\ \mathrm{(g\,cm^{-3})}$')
ax.set_ylabel(r'$\log_{10}P\ \mathrm{(dyne\,cm^{-2})}$')
ax.set_title(f"Composite EOS (SBD + BBP + Magnetic BPS, b={B_field_str}) with C¹ joins")
ax.legend()
ax.grid(True, alpha=0.3)
plt.show()
