
import math, re
from pathlib import Path

ROOT = Path(__file__).resolve().parent
BASE = ROOT / "couette-3d.prm"
EPSILONS = (0.002, 0.005, 0.008)

src = BASE.read_text(encoding="utf-8")
for eps in EPSILONS:
    radius = (eps * 75 * 30 * 20 / (10 * 4 / 3 * math.pi)) ** (1 / 3)
    tag = f"eps{round(eps * 1000):03d}"
    text = re.sub(r"(set shape arguments\s*=\s*)([0-9.eE+-]+)", rf"\g<1>{radius:.8f}", src)
    text = re.sub(r"(set output name\s*=\s*)\S+", rf"\g<1>couette_3d_{tag}", text)
    text = re.sub(r"(set force name\s*=\s*)\S+", rf"\g<1>force_3d_{tag}", text)
    text = re.sub(r"(set output path\s*=\s*)\S+", rf"\g<1>out_{tag}/", text)
    folder = ROOT / tag
    folder.mkdir(parents=True, exist_ok=True)
    (folder / f"couette_3d_{tag}.prm").write_text(text, encoding="utf-8")
    print(f"{eps:.1%}   rayon = {radius:.4f}   Einstein = {10 * (1 + 2.5 * eps):.3f}")
