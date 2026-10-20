# SPDX-FileCopyrightText: Copyright (c) 2023-2026 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
"""Generate 2D Couette .prm files for several solid fractions phi.

Reads couette-2d-corr.prm and, for each phi in EPSILONS, writes a new
parameter file with the matching particle radius and output names.
"""

import math
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent
BASE = ROOT / "couette-2d-corr.prm"  # base 2D Couette parameter file
LX, LY, N = 50.0, 20.0, 7            # domain size and number of discs
EPSILONS = (0.02, 0.04, 0.05, 0.15)  # solid fractions phi to generate

src = BASE.read_text(encoding="utf-8")
for eps in EPSILONS:
    # radius from phi = N * pi * R^2 / (LX * LY)
    radius = math.sqrt(eps * LX * LY / (N * math.pi))
    tag = f"eps{round(eps * 100):02d}"

    # copy base prm and update radius + output/force/path names
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
    print(f"{eps:.0%}   rayon = {radius:.4f}   {folder.name}/{prm.name}")
