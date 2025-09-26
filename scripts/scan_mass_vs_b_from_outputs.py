import re, csv, argparse
from pathlib import Path
from math import isfinite

M_SUN = 1.98847e33  # g

# b directory: b_0e00, b_1e-05, b_1e-03, ...
B_DIR_RX = re.compile(r"^b_([0-9.+\-eE]+)$")

# central density: trailing scientific number before extension
# ..._1.00e+14.csv  -> capture "1.00e+14"
TRAILING_SCI_RX = re.compile(r"_([0-9]+(?:\.[0-9]+)?[eE][+\-]?[0-9]+)\.(csv|dat)$")

def parse_args():
    ap = argparse.ArgumentParser(
        description="Build M(b) at fixed central density using TOV outputs in data/b_* folders."
    )
    ap.add_argument("--base", default="data", help="base directory containing b_* folders")
    ap.add_argument("--rho_c", type=float, required=True, help="target central density [g/cm^3]")
    ap.add_argument("--out", default="data/mass_vs_b.csv", help="output CSV file")
    ap.add_argument("--ext", default=".csv", help="file extension to scan")
    ap.add_argument("--strict-common", action="store_true",
                    help="require a common rho_c across all b folders (use the one closest to --rho_c)")
    ap.add_argument("--rtol", type=float, default=5e-4,
                    help="relative tolerance to match common rho_c across folders (default=5e-4)")
    return ap.parse_args()

def list_b_folders(base: Path):
    for d in sorted(p for p in base.iterdir() if p.is_dir()):
        m = B_DIR_RX.match(d.name)
        if m:
            yield float(m.group(1)), d

def parse_rho_from_name(p: Path):
    m = TRAILING_SCI_RX.search(p.name)
    if not m:
        return None
    try:
        return float(m.group(1))
    except ValueError:
        return None

def gather_folder_rhos(folder: Path, ext: str):
    vals = []
    for f in folder.glob(f"*{ext}"):
        rho = parse_rho_from_name(f)
        if rho is not None and isfinite(rho):
            vals.append((rho, f))
    return sorted(vals, key=lambda t: t[0])

def pick_nearest(target: float, pairs):
    """pairs: list of (rho, Path)"""
    best = None
    for rho, f in pairs:
        d = abs(rho - target)
        if best is None or d < best[0]:
            best = (d, rho, f)
    return best  # (delta, rho, path) or None

def group_rhos(values, rtol):
    """
    values: iterable of floats
    returns list of groups, each group is list of floats that agree within rtol
    """
    vals = sorted(values)
    groups = []
    cur = []
    for v in vals:
        if not cur:
            cur = [v]
        else:
            if abs(v - cur[-1]) <= rtol * max(abs(v), abs(cur[-1]), 1.0):
                cur.append(v)
            else:
                groups.append(cur)
                cur = [v]
    if cur:
        groups.append(cur)
    return groups

def pick_common_rho(all_sets, target, rtol):
    """
    all_sets: list of sorted unique rho arrays per folder
    returns rho_common or None
    """
    # merge and group all rhos
    union = sorted({rho for s in all_sets for rho in s})
    groups = group_rhos(union, rtol)
    # keep only groups present in every folder (by rtol match)
    candidates = []
    for g in groups:
        rep = sum(g) / len(g)  # representative rho
        ok = True
        for s in all_sets:
            # any s member within rtol of rep?
            if not any(abs(r - rep) <= rtol * max(abs(r), abs(rep), 1.0) for r in s):
                ok = False
                break
        if ok:
            candidates.append(rep)
    if not candidates:
        return None
    # choose candidate closest to target
    return min(candidates, key=lambda r: abs(r - target))

def mass_last_row_Msun(csv_path: Path) -> float:
    # Format: log_r[cm],log_m[g],log_P[dyne/cm^2]
    last = None
    with csv_path.open() as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            last = s
    if last is None or "log_" in last:
        raise RuntimeError(f"No data rows in {csv_path}")
    cols = [c.strip() for c in last.split(",")]
    if len(cols) < 2:
        raise RuntimeError(f"Unexpected row in {csv_path}: {last}")
    log_m_g = float(cols[1])
    return (10.0 ** log_m_g) / M_SUN

def main():
    args = parse_args()
    base = Path(args.base)

    # scan each b folder
    per_b = []  # [(b, [(rho, path), ...]), ...]
    for b_val, folder in list_b_folders(base):
        pairs = gather_folder_rhos(folder, args.ext)
        if not pairs:
            print(f"[skip] no files with trailing sci-number in {folder}")
            continue
        per_b.append((b_val, pairs))

    if not per_b:
        print("No data found.")
        return

    rows = []
    if args.strict_common:
        # compute common rho across all b's (within rtol), pick nearest to target
        sets = [sorted({rho for rho, _ in pairs}) for _, pairs in per_b]
        rho_common = pick_common_rho(sets, args.rho_c, args.rtol)
        if rho_common is None:
            print("[strict-common] no common ρc across all b folders within tolerance.")
            return
        print(f"[strict-common] using common ρc ≈ {rho_common:.6e}")
        for b_val, pairs in per_b:
            _, rho_used, fpath = pick_nearest(rho_common, pairs)
            M = mass_last_row_Msun(fpath)
            rows.append({"b": b_val, "rho_c_target": rho_common,
                        "rho_c_used": rho_used, "delta_rho": abs(rho_used - rho_common),
                        "M_Msun": M, "file": str(fpath)})
            print(f"[ok] b={b_val:g} -> {fpath.name} (ρc≈{rho_used:.3e}) M={M:.5f} Msun")
    else:
        # per-b nearest to user target
        for b_val, pairs in per_b:
            pick = pick_nearest(args.rho_c, pairs)
            if pick is None:  # defensive
                print(f"[skip] no parseable ρc in {b_val}")
                continue
            delta, rho_used, fpath = pick
            M = mass_last_row_Msun(fpath)
            rows.append({"b": b_val, "rho_c_target": args.rho_c,
                        "rho_c_used": rho_used, "delta_rho": delta,
                        "M_Msun": M, "file": str(fpath)})
            print(f"[ok] b={b_val:g} -> {fpath.name} (ρc≈{rho_used:.3e}) M={M:.5f} Msun")

    # sort and write
    rows.sort(key=lambda r: r["b"])
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["b","rho_c_target","rho_c_used","delta_rho","M_Msun","file"])
        w.writeheader()
        w.writerows(rows)

    print("Wrote", out)

if __name__ == "__main__":
    main()
