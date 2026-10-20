# SPDX-FileCopyrightText: Copyright (c) 2023-2026 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
"""Generate 3D Couette .prm files for several solid fractions phi.

Reads couette-3d-20p.prm and, for each phi in EPSILONS, writes a new
parameter file with the matching particle radius and output names.
"""

import math
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent
BASE = ROOT / "couette-3d-20p.prm"       # base 3D Couette parameter file
LX, LY, LZ, N = 75.0, 30.0, 20.0, 20      # domain size and number of spheres
EPSILONS = (0.002, 0.005, 0.008)         # solid fractions phi to generate

src = BASE.read_text(encoding="utf-8")
for eps in EPSILONS:
    # radius from phi = N * (4/3) * pi * R^3 / (LX * LY * LZ)
    radius = (eps * LX * LY * LZ / (N * 4.0 / 3.0 * math.pi)) ** (1.0 / 3.0)
    tag = f"eps{round(eps * 1000):03d}"

    # copy base prm and update radius + output/force/path names
    text = re.sub(
        r"(set shape arguments\s*=\s*)([0-9.eE+-]+)",
        rf"\g<1>{radius:.8f}",
        src,
    )
    text = re.sub(r"(set output name\s*=\s*)\S+", rf"\g<1>couette_3d_20p_{tag}", text)
    text = re.sub(r"(set force name\s*=\s*)\S+", rf"\g<1>force_20p_{tag}", text)
    text = re.sub(r"(set output path\s*=\s*)\S+", rf"\g<1>out_20p_{tag}/", text)

    folder = ROOT / f"20p_{tag}"
    folder.mkdir(parents=True, exist_ok=True)
    prm = folder / f"couette_3d_20p_{tag}.prm"
    prm.write_text(text, encoding="utf-8")
    print(f"{eps:.1%}   R={radius:.4f}   {prm}")
