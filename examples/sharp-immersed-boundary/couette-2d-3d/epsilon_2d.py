# SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
"""Generate 2D Couette .prm files for several solid fractions phi.

Reads couette-2d-7p.prm and, for each phi in EPSILONS, writes a new
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
# Channel 50 x 20 with 7 discs. Only the radius changes with phi.
ROOT = Path(__file__).resolve().parent
BASE = ROOT / "couette-2d-7p.prm"
LX = 50.0
LY = 20.0
N = 7
EPSILONS = (0.02, 0.04, 0.05, 0.15)

#---------------------------------------------------
# Write One Parameter File per Solid Fraction
#---------------------------------------------------
# Area fraction:  phi = N * pi * R^2 / (LX * LY)
# Radius:         R   = sqrt(phi * LX * LY / (N * pi))
src = BASE.read_text(encoding="utf-8")

for eps in EPSILONS:
    radius = math.sqrt(eps * LX * LY / (N * math.pi))
    tag = f"eps{round(eps * 100):02d}"

    # Copy the base prm and replace radius, output name, force name, path.
    text = re.sub(
        r"(set shape arguments\s*=\s*)([0-9.eE+-]+)",
        rf"\g<1>{radius:.8f}",
        src,
    )
    text = re.sub(r"(set output name\s*=\s*)\S+", rf"\g<1>couette_{tag}", text)
    text = re.sub(r"(set force name\s*=\s*)\S+", rf"\g<1>force_{tag}", text)
    text = re.sub(r"(set output path\s*=\s*)\S+", rf"\g<1>out_{tag}/", text)

    folder = ROOT / tag
    folder.mkdir(parents=True, exist_ok=True)
    prm = folder / f"couette_{tag}.prm"
    prm.write_text(text, encoding="utf-8")
    print(f"{eps:.0%}   radius = {radius:.4f}   {folder.name}/{prm.name}")
