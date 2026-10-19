
import math
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent
BASE = ROOT / "couette-3d-20p.prm"
LX, LY, LZ, N = 75.0, 30.0, 20.0, 20
EPSILONS = (0.002, 0.005, 0.008)

src = BASE.read_text(encoding="utf-8")
for eps in EPSILONS:
    radius = (eps * LX * LY * LZ / (N * 4.0 / 3.0 * math.pi)) ** (1.0 / 3.0)
    tag = f"eps{round(eps * 1000):03d}"
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
