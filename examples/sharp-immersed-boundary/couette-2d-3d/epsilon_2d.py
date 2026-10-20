# SPDX-FileCopyrightText: Copyright (c) 2023-2026 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later


import math, re
from pathlib import Path

ROOT = Path(__file__).resolve().parent
BASE = ROOT / "couette-2d-corr.prm"
LX, LY, N = 50.0, 20.0, 7
EPSILONS = (0.02, 0.04, 0.05, 0.15)

src = BASE.read_text(encoding="utf-8")
for eps in EPSILONS:
    radius = math.sqrt(eps * LX * LY / (N * math.pi))
    tag = f"eps{round(eps * 100):02d}"
    text = re.sub(r"(set shape arguments\s*=\s*)([0-9.eE+-]+)", rf"\g<1>{radius:.8f}", src)
    text = re.sub(r"(set output name\s*=\s*)\S+", rf"\g<1>couette_{tag}", text)
    text = re.sub(r"(set force name\s*=\s*)\S+", rf"\g<1>force_{tag}", text)
    text = re.sub(r"(set output path\s*=\s*)\S+", rf"\g<1>out_{tag}/", text)
    folder = ROOT / tag
    folder.mkdir(parents=True, exist_ok=True)
    prm = folder / f"couette_{tag}.prm"
    prm.write_text(text, encoding="utf-8")
    print(f"{eps:.0%}   rayon = {radius:.4f}   {folder.name}/{prm.name}")
