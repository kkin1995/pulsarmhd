#!/usr/bin/env python3
import os, glob, re, argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

SOLAR_MASS_G = 1.989e33
CM_TO_KM     = 1e-5

def list_b_dirs(root):
    """Return [(b_value, path), ...] sorted by b."""
    pairs = []
    for d in sorted(glob.glob(os.path.join(root, "b_*"))):
        base = os.path.basename(d)           # e.g. b_1e-03
        m = re.match(r"b_([0-9.eE+-]+)$", base)
        if not m: 
            continue
        try:
            b = float(m.group(1))
        except ValueError:
            continue
        pairs.append((b, d))
    return sorted(pairs, key=lambda x: x[0])

def read_surface_mr(csv_path, r_end_cm=None):
    """Return (M_solar, R_km) from last row, or None if invalid."""
    df = pd.read_csv(csv_path)
    needed = {'log_m[g]', 'log_r[cm]'}
    if not needed.issubset(df.columns):
        return None

    # optional quick surface check if pressure column present
    if 'log_P[dyne/cm^2]' in df.columns and len(df) > 3:
        lp = df['log_P[dyne/cm^2]'].to_numpy()
        if not np.all(np.diff(lp) <= 1e-6):   # allow tiny wiggles
            return None
        # last pressure ~ minimum pressure
        if not np.isclose(lp[-1], lp.min(), atol=5e-4):
            return None

    m_sol = 10.0**df['log_m[g]'].iloc[-1] / SOLAR_MASS_G
    r_km  = 10.0**df['log_r[cm]'].iloc[-1] * CM_TO_KM

    # skip models that likely hit outer boundary instead of the surface
    if r_end_cm is not None:
        r_end_km = r_end_cm * CM_TO_KM
        if (r_km > 0.98*r_end_km) and np.isclose(r_km, r_end_km, rtol=1e-3, atol=0.5):
            return None

    return m_sol, r_km

def collect_family(b_dir, r_end_cm=None):
    """Return arrays (R_km_sorted, M_solar_sorted) from all CSVs in b_dir."""
    all_pairs = []
    for csv in glob.glob(os.path.join(b_dir, "*.csv")):
        mr = read_surface_mr(csv, r_end_cm=r_end_cm)
        if mr is not None:
            all_pairs.append(mr)
    if not all_pairs:
        return np.array([]), np.array([])
    M, R = zip(*all_pairs)
    M = np.array(M); R = np.array(R)
    # sort by radius for a tidy curve
    idx = np.argsort(R)
    return R[idx], M[idx]

def main():
    ap = argparse.ArgumentParser(description="Quick M(R) plot colored by b (reads data/b_*/*csv)")
    ap.add_argument("--data-dir", default="data", help="top-level data directory (default: data)")
    ap.add_argument("--out", default="mr_by_b.png", help="output figure path")
    ap.add_argument("--r-end-cm", type=float, default=1.0e8,
                    help="if final R ~ r_end, treat as non-surface and skip (default 1e8 cm)")
    args = ap.parse_args()

    b_dirs = list_b_dirs(args.data_dir)
    if not b_dirs:
        print("No data/b_* directories found.")
        return

    ref_b = None                 # set after we find the smallest b
    ref_curve = None             # (R_km, M_solar)

    fig = plt.figure(figsize=(7.5, 7.0), constrained_layout=True)
    gs = fig.add_gridspec(2, 1, height_ratios=[3, 1])
    ax  = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1], sharex=ax)

    summary = []
    for b, path in b_dirs:
        files = sorted(glob.glob(os.path.join(path, "*.csv")))
        R_km, M_solar = collect_family(path, r_end_cm=args.r_end_cm)
        used = len(R_km)
        summary.append((b, len(files), used, 
                        (float(np.min(R_km)) if used else np.nan,
                        float(np.max(R_km)) if used else np.nan)))

        if used == 0:
            print(f"[warn] {path}: {len(files)} files, 0 usable (likely r_end/surface filter)")
            continue

        label = f"b = {b:g}"
        #ax.plot(R_km, M_solar, marker="o", ms=2.6, lw=1.3, label=label, alpha=0.9)
        ax.plot(R_km, M_solar, marker="o", ms=2.6, lw=1.6, alpha=0.85, zorder=2, label=label)
        ax.grid(True, alpha=0.3, zorder=0)


        # pick reference as the smallest b we actually plotted
        if ref_curve is None or (ref_b is not None and b < ref_b):
            ref_b, ref_curve = b, (R_km, M_solar)

    # residuals vs reference
    if ref_curve is not None:
        Rref, Mref = ref_curve
        # build an interpolator in M(R) for the reference
        from scipy.interpolate import interp1d
        fMref = interp1d(Rref, Mref, bounds_error=False, fill_value=np.nan)

        # overlay residuals for each b (skip the reference itself)
        for b, path in b_dirs:
            R_km, M_solar = collect_family(path, r_end_cm=args.r_end_cm)
            if R_km.size == 0 or b == ref_b: 
                continue
            # evaluate ΔM at radii common to both (clip to overlap range)
            Rmin = max(np.min(R_km), np.min(Rref))
            Rmax = min(np.max(R_km), np.max(Rref))
            mask = (R_km >= Rmin) & (R_km <= Rmax)
            dM = M_solar[mask] - fMref(R_km[mask])
            ax2.plot(R_km[mask], dM, lw=1.2, label=f"ΔM (b={b:g} − {ref_b:g})")

    # cosmetics
    ax.set_title(r"Mass–Radius curves colored by $b = B/B_{\rm Q}$")
    ax.set_xlabel("Radius R (km)")
    ax.set_ylabel(r"Gravitational Mass M ($M_\odot$)")
    ax.grid(True, alpha=0.3)
    ax.legend(frameon=True, fontsize=9)

    ax2.axhline(0, ls=":", lw=1)
    ax2.set_xlabel("Radius R (km)")
    ax2.set_ylabel(r"ΔM ($M_\odot$)") 
    ax2.grid(True, alpha=0.3)

    plt.savefig(args.out, dpi=200)
    print(f"Saved {args.out}")

    # log a compact per-b summary
    print("\nPer-b summary (files, used, R-range [km]):")
    for b, nfiles, used, Rrng in summary:
        rmin, rmax = Rrng
        print(f"  b={b:g}: files={nfiles:3d}, used={used:3d}, "
            f"R∈[{'' if np.isnan(rmin) else f'{rmin:.1f}'}, "
            f"{'' if np.isnan(rmax) else f'{rmax:.1f}'}]")


if __name__ == "__main__":
    main()
