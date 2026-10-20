
# SPDX-FileCopyrightText: Copyright (c) 2023-2026 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
"""Generate 20-particle center arrangements D1–D4 for the 3D Couette example.

Domain size matches couette-3d-20p.prm. Centers lie in the mid-plane z = Z.
Each .dat line is: x y z vx, with vx = y / LY (local Couette velocity).
"""

import math
import random
from pathlib import Path

# LX, LY, LZ: domain length in x, y, and z
LX, LY, LZ = 75.0, 30.0, 20.0
# Z: z-coordinate of all particle centers (mid-plane)
# R: particle radius
# N: number of particles
Z, R, N = 10.0, 1.0, 20
# OUT: output directory for arrangement files
OUT = Path("arrangements")
OUT.mkdir(exist_ok=True)


def write(name, pts):
    # Write particle centers to arrangements/<name>.dat
    # Each line: x y z vx  (vx = Couette velocity at the center)
    lines = [f"# {name}  N={len(pts)}  z={Z}"]
    for x, y, z in pts:
        lines.append(f"{x:.8f} {y:.8f} {z:.8f} {y / LY:.8f}")
    path = OUT / f"{name}.dat"
    path.write_text("\n".join(lines) + "\n")
    print(name, len(pts), "->", path)


def D1():
    # Regular grid: 5 columns x 4 rows
    xs = [7.5 + i * 15.0 for i in range(5)]
    ys = [4.5, 11.5, 18.5, 25.5]
    return [(x, y, Z) for y in ys for x in xs]


def D2():
    # Two rows: 10 particles at y=8 and 10 at y=22
    xs = [5.0 + i * (65.0 / 9.0) for i in range(10)]
    return [(x, 8.0, Z) for x in xs] + [(x, 22.0, Z) for x in xs]


def D3():
    # Staggered grid: like D1, odd rows shifted in x
    xs_a = [7.5 + i * 15.0 for i in range(5)]
    xs_b = [15.0 + i * 15.0 for i in range(5)]
    pts = []
    for j, y in enumerate([4.5, 11.5, 18.5, 25.5]):
        xs = xs_a if j % 2 == 0 else xs_b
        pts.extend((x, y, Z) for x in xs)
    return pts


def D4():
    # Random positions with fixed seed and minimum gap between centers
    random.seed(7)
    pts, gap = [], 5.0
    while len(pts) < N:
        p = (random.uniform(3, 72), random.uniform(3, 27), Z)
        if all(math.hypot(p[0] - q[0], p[1] - q[1]) >= gap for q in pts):
            pts.append(p)
    return pts


for name, fn in [("D1", D1), ("D2", D2), ("D3", D3), ("D4", D4)]:
    write(name, fn())
