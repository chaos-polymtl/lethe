# SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
"""Generate 3D Couette .prm files for several solid fractions phi.

Reads couette-3d-20p.prm and, for each phi in EPSILONS, writes a new
parameter file with the matching particle radius and output names.
"""

#---------------------------------------------------
# Imports
#---------------------------------------------------
import math
import re
from pathlib import Path

#---------------------------------------------------
# Base File and Geometry
#---------------------------------------------------
# Box 75 x 30 x 20 with 20 spheres. Only the radius changes with phi.
# The same midplane centers overlap above about 1 %, so the sweep
# stays at 0.2 %, 0.5 % and 0.8 %.
ROOT = Path(__file__).resolve().parent
BASE = ROOT / "couette-3d-20p.prm"
LX = 75.0
LY = 30.0
LZ = 20.0
N = 20
EPSILONS = (0.002, 0.005, 0.008)

#---------------------------------------------------
# Write One Parameter File per Solid Fraction
#---------------------------------------------------
# Volume fraction:  phi = N * (4/3) * pi * R^3 / (LX * LY * LZ)
# Radius:           R   = [phi * LX * LY * LZ / (N * 4/3 * pi)]^(1/3)
src = BASE.read_text(encoding="utf-8")

for eps in EPSILONS:
    radius = (eps * LX * LY * LZ / (N * 4.0 / 3.0 * math.pi)) ** (1.0 / 3.0)
    tag = f"eps{round(eps * 1000):03d}"

    # Copy the base prm and replace radius, output name, force name, path.
    text = re.sub(
        r"(set shape arguments\s*=\s*)([0-9.eE+-]+)",
        rf"\g<1>{radius:.8f}",
        src,
    )
    text = re.sub(
        r"(set output name\s*=\s*)\S+",
        rf"\g<1>couette_3d_20p_{tag}",
        text,
    )
    text = re.sub(r"(set force name\s*=\s*)\S+", rf"\g<1>force_20p_{tag}", text)
    text = re.sub(r"(set output path\s*=\s*)\S+", rf"\g<1>out_20p_{tag}/", text)

    folder = ROOT / f"20p_{tag}"
    folder.mkdir(parents=True, exist_ok=True)
    prm = folder / f"couette_3d_20p_{tag}.prm"
    prm.write_text(text, encoding="utf-8")
    print(f"{eps:.1%}   radius = {radius:.4f}   {folder.name}/{prm.name}")
