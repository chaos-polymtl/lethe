import math
import random
from pathlib import Path

LX, LY, LZ = 75.0, 30.0, 20.0
Z, R, N = 10.0, 1.0, 20
OUT = Path("arrangements")
OUT.mkdir(exist_ok=True)

def write(name, pts):
    lines = [f"# {name}  N={len(pts)}  z={Z}"]
    for x, y, z in pts:
        lines.append(f"{x:.8f} {y:.8f} {z:.8f} {y/20.0:.8f}")
    path = OUT / f"{name}.dat"
    path.write_text("\n".join(lines) + "\n")
    print(name, len(pts), "->", path)

def D1():
    xs = [7.5 + i * 15.0 for i in range(5)]
    ys = [4.5, 11.5, 18.5, 25.5]
    return [(x, y, Z) for y in ys for x in xs]

def D2():
    xs = [5.0 + i * (65.0 / 9.0) for i in range(10)]
    return [(x, 8.0, Z) for x in xs] + [(x, 22.0, Z) for x in xs]

def D3():
    xs_a = [7.5 + i * 15.0 for i in range(5)]
    xs_b = [15.0 + i * 15.0 for i in range(5)]
    pts = []
    for j, y in enumerate([4.5, 11.5, 18.5, 25.5]):
        xs = xs_a if j % 2 == 0 else xs_b
        pts.extend((x, y, Z) for x in xs)
    return pts

def D4():
    random.seed(7)
    pts, gap = [], 5.0
    while len(pts) < N:
        p = (random.uniform(3, 72), random.uniform(3, 27), Z)
        if all(math.hypot(p[0] - q[0], p[1] - q[1]) >= gap for q in pts):
            pts.append(p)
    return pts

for name, fn in [("D1", D1), ("D2", D2), ("D3", D3), ("D4", D4)]:
    write(name, fn())
