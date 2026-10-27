# SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
"""Plot 2D Couette mu_eff versus the dilute disc law mu (1 + 2 phi).

The empty-channel top-wall force is F0 = 25. With particles the last
recorded |Fx| gives mu_eff = MU0 * |Fx| / F0. The reference curve uses
[eta] = 2 (circular cylinders), not the 3D Einstein coefficient 5/2.
"""

#---------------------------------------------------
# Imports
#---------------------------------------------------
from pathlib import Path
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
DOC = Path.home() / "lethe/doc/source/examples/sharp-immersed-boundary/couette-2d-3d/images"

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 12,
    "axes.labelsize": 14,
    "legend.fontsize": 11,
})

#---------------------------------------------------
# Physical Parameters
#---------------------------------------------------
# mu = 10, U = 1, h = 20, Lx = 50
# F0 = mu * U * Lx / h = 25
MU0 = 10.0
F0 = 25.0
EPSILONS = (0.02, 0.04, 0.05, 0.15)


def get_fx(filename):
    """Return |f_x| from the last row of a Lethe force file."""
    with open(filename, "r") as f:
        lines = [line.strip() for line in f if line.strip()]
    col = lines[0].split().index("f_x")
    return abs(float(lines[-1].split()[col]))


def find_force_file(tag):
    """Locate the top-wall force file written for one solid fraction."""
    candidates = [
        ROOT / f"out_{tag}" / f"force_{tag}.03.dat",
        ROOT / tag / f"out_{tag}" / f"force_{tag}.03.dat",
        ROOT / tag / f"force_{tag}.03.dat",
    ]
    for path in candidates:
        if path.exists():
            return path
    return None


#---------------------------------------------------
# Effective Viscosity from Simulation
#---------------------------------------------------
phi_pct, mu_eff = [], []
for eps in EPSILONS:
    tag = f"eps{round(eps * 100):02d}"
    path = find_force_file(tag)
    if path is None:
        print(f"MISSING phi={100 * eps:.0f}%")
        continue
    Fx = get_fx(path)
    mu = MU0 * Fx / F0
    ein = MU0 * (1.0 + 2.0 * eps)
    phi_pct.append(100 * eps)
    mu_eff.append(mu)
    print(f"phi={100 * eps:.0f}%  |Fx|={Fx:.4f}  mu_eff={mu:.4f}  Einstein={ein:.4f}")

if len(phi_pct) < 2:
    raise SystemExit("Not enough force files")

#---------------------------------------------------
# Dilute 2D Reference  mu (1 + 2 phi)
#---------------------------------------------------
phi_max = max(phi_pct) * 1.05
n = 200
phi_line = [phi_max * i / (n - 1) for i in range(n)]
mu_line = []
for p in phi_line:
    mu_line.append(MU0 * (1.0 + 2.0 * p / 100.0))

#---------------------------------------------------
# Figure
#---------------------------------------------------
plt.figure(figsize=(6, 4))
plt.plot(phi_line, mu_line, "k-", label=r"Einstein $\mu\left(1+2\phi\right)$")
plt.plot(phi_pct, mu_eff, "o", color="#1f4e79", label=r"Simulation $\mu_{\mathrm{eff}}$")
plt.xlabel(r"Solid fraction $\phi$ (%)")
plt.ylabel(r"Effective viscosity $\mu_{\mathrm{eff}}$")
plt.legend()
plt.xlim(left=0.0)
plt.tight_layout()
DOC.mkdir(parents=True, exist_ok=True)
plt.savefig(ROOT / "viscosity-2d.png", dpi=150)
plt.savefig(DOC / "viscosity-2d.png", dpi=150)
print("wrote viscosity-2d.png")
