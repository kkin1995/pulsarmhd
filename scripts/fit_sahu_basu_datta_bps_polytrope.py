import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator
import matplotlib.pyplot as plt
import math

speed_of_light = 29979245800  # cm s-1

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

# ---------- Load EOS pieces ----------
sahu_basu_datta_eos = pd.read_csv(
    "/home/karan-kinariwala/Dropbox/KARAN/2-Areas/Education/PhD/3-Research/pulsarmhd/data/sahu-basu-datta.dat",
    sep=r'\s+', header=None, dtype=float
)
print(f"Number of NaN Values in Sahu Basu Datta: {sahu_basu_datta_eos.isna().sum(axis = 0)}")

outfile = "data/sahu_basu_datta_bbp_magnetic_bps_b_1e-05.csv"

magnetic_bps_eos = pd.read_csv(
    "/home/karan-kinariwala/Dropbox/KARAN/2-Areas/Education/PhD/3-Research/pulsarmhd/data/bps_single_b_b_1e-05.csv",
    header=0
)
print(f"Number of NaN Values in Magnetic BPS: {magnetic_bps_eos.isna().sum(axis = 0)}")

bbp_eos = pd.read_csv(
    "/home/karan-kinariwala/Dropbox/KARAN/2-Areas/Education/PhD/3-Research/pulsarmhd/data/bbp.csv",
    header = 0
)
print(f"Number of NaN Values in BBP: {bbp_eos.isna().sum(axis = 0)}")

# ---------- Splicing windows ----------
neutron_drip = 4.46e11
sbd_mass_density_min = sahu_basu_datta_eos.iloc[:, 1].min() * 1e14
bbp_mass_density_max = bbp_eos.iloc[:, 0].max()

# BPS: replace ±inf with NaN and drop, then keep up to drip (log-space)
for col in ['log_n', 'log_rho', 'log_P']:
    magnetic_bps_eos[col] = magnetic_bps_eos[col].replace([np.inf, -np.inf], np.nan)
magnetic_bps_eos = magnetic_bps_eos.dropna(subset=['log_rho','log_P'])
magnetic_bps_eos = magnetic_bps_eos[magnetic_bps_eos['log_rho'].to_numpy() <= np.log10(neutron_drip)]

# BBP: keep in [drip, bbp_max]
bbp_eos = bbp_eos[(bbp_eos.iloc[:, 0] >= neutron_drip) & (bbp_eos.iloc[:, 0] <= bbp_mass_density_max)]

# SBD: keep above BBP max (its table is in units of 1e14)
sahu_basu_datta_eos = sahu_basu_datta_eos[sahu_basu_datta_eos.iloc[:, 1] >= (bbp_mass_density_max / 1e14)]

# ---------- Transform to logs ----------
sbd_number_density = sahu_basu_datta_eos.iloc[:, 0].values * 1e14
sbd_mass_density   = sahu_basu_datta_eos.iloc[:, 1].values * 1e14
sbd_pressure       = sahu_basu_datta_eos.iloc[:, 2].values * (speed_of_light ** 2) * 1e14
print(f"Density Range (Sahu Basu Datta): Min: {min(sbd_mass_density)} g cm-3 | Max: {max(sbd_mass_density)} g cm-3")

sbd_log_n   = np.log10(sbd_number_density)
sbd_log_rho = np.log10(sbd_mass_density)
sbd_log_P   = np.log10(sbd_pressure)

magnetic_bps_log_rho = magnetic_bps_eos.iloc[:, 1].values
magnetic_bps_log_P   = magnetic_bps_eos.iloc[:, 2].values
print(f"Density Range (Magnetic BPS): Min: {min(10 ** magnetic_bps_log_rho)} g cm-3 | Max: {max(10 ** magnetic_bps_log_rho)} g cm-3")

# BBP to logs (guard positivity)
bbp_mass_density = bbp_eos.iloc[:, 0].values
bbp_pressure     = bbp_eos.iloc[:, 1].values
mask_bbp = np.isfinite(bbp_mass_density) & np.isfinite(bbp_pressure) & (bbp_mass_density > 0) & (bbp_pressure > 0)
bbp_mass_density = bbp_mass_density[mask_bbp]
bbp_pressure     = bbp_pressure[mask_bbp]
bbp_log_rho = np.log10(bbp_mass_density)
bbp_log_P   = np.log10(bbp_pressure)
print(f"Density Range (BBP): Min: {min(bbp_mass_density)} g cm-3 | Max: {max(bbp_mass_density)} g cm-3")

# ---------- Minimal cleaning for PCHIP ----------
magnetic_bps_log_rho, magnetic_bps_log_P = clean_monotone_xy(magnetic_bps_log_rho, magnetic_bps_log_P)
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

pd.DataFrame({"log_rho": comp_log_rho, "log_P": comp_log_P}).to_csv(outfile, index=False)

# ---------- Diagnostics ----------
print(f"Magnetic BPS max: log(ρ)={magnetic_bps_log_rho[-1]:.2f}, log(P)={magnetic_bps_log_P[-1]:.2f}")
print(f"BBP min:         log(ρ)={bbp_log_rho[0]:.2f},  log(P)={bbp_log_P[0]:.2f}")
print(f"BBP max:         log(ρ)={bbp_log_rho[-1]:.2f}, log(P)={bbp_log_P[-1]:.2f}")
print(f"SBD min:         log(ρ)={sbd_log_rho[0]:.2f},  log(P)={sbd_log_P[0]:.2f}")

# ---------- Plot ----------
fig, ax = plt.subplots(figsize=(10, 5))
ax.scatter(sbd_log_rho, sbd_log_P, s=5, c='red',    label="Sahu–Basu–Datta (core, trimmed)")
ax.scatter(bbp_log_rho, bbp_log_P, s=5, c='green',  label="Baym–Bethe–Pethick (inner crust, trimmed)")
ax.scatter(magnetic_bps_log_rho, magnetic_bps_log_P, s=5, c='purple', label="Magnetic BPS (outer crust, b=1e-1)")
ax.plot(xb1, yb1, 'k.', ms=6, label='C¹ bridge (BPS→BBP)')
ax.plot(xb2, yb2, 'k.', ms=6, label='C¹ bridge (BBP→SBD)')
ax.set_xlabel(r'$\log_{10}\rho\ \mathrm{(g\,cm^{-3})}$')
ax.set_ylabel(r'$\log_{10}P\ \mathrm{(dyne\,cm^{-2})}$')
ax.set_title("Composite EOS (SBD + BBP + Magnetic BPS, b=1e-1) with C¹ joins")
ax.legend()
ax.grid(True, alpha=0.3)
plt.show()
