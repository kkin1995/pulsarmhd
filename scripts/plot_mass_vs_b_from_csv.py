# plot_mass_vs_b_from_csv
import csv
import matplotlib.pyplot as plt
from pathlib import Path

in_csv = Path("data/mass_vs_b.csv")

# --- 1) Read b exactly as stored (including 0.0) ---
b_vals, M_vals = [], []
with in_csv.open() as f:
    rd = csv.DictReader(f)
    for r in rd:
        b_vals.append(float(r["b"]))       # do NOT replace 0 with epsilon
        M_vals.append(float(r["M_Msun"]))

fig, ax = plt.subplots()

# --- 2) Use symlog so x=0 is valid; keep a tiny linear window around 0 ---
linthresh = 1e-12
ax.set_xscale("symlog", linthresh=linthresh)

# Plot
ax.plot(b_vals, M_vals, marker="o")

# --- 3) Explicit ticks so 0 is shown as '0' and no phantom 1e-11 tick appears ---
pos_bs = sorted({b for b in b_vals if b > 0})
ticks   = [0.0] + pos_bs
labels  = ["0"] + [f"{b:g}" for b in pos_bs]
ax.set_xticks(ticks)
ax.set_xticklabels(labels)

# Limits: small padding; left bound slightly negative so the 0 marker isn’t clipped
right = (max(pos_bs) if pos_bs else linthresh) * 1.2
ax.set_xlim(-linthresh, right)

ax.set_xlabel(r"$b = B / B_Q$")
ax.set_ylabel(r"Gravitational mass $M$ ($M_\odot$)")
ax.grid(True, which="both", linestyle=":")

fig.tight_layout()
plt.show()
