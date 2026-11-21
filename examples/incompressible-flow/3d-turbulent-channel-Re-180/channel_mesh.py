#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

"""Mesh of the turbulent channel: the grid arguments of channel.prm.

The coarse mesh is deal.II's subdivided_hyper_rectangle with step sizes:
uniform cells in x and z, and cells stretched towards the walls in y,

    y_j = tanh(gamma (2 j / n_y - 1)) / tanh(gamma),   j = 0 .. n_y,

walls at y = -1 and 1 (delta = 1). Each initial refinement halves every cell,
so the refined meshes keep the stretching and are nested.

Prints the grid arguments and the wall resolution of each refinement, in the
wall units of the DNS (Re_tau = 182). With --prm, also writes the grid
arguments into that parameter file.

The defaults are the mesh of channel.prm: a 4 pi x 2 x 2 pi box with
32 x 16 x 32 coarse cells and gamma = 2.

Usage:
  python3 channel_mesh.py                          # the mesh of channel.prm
  python3 channel_mesh.py --gamma 2.2 --prm channel.prm
  python3 channel_mesh.py --cells 16 16 16 --lx 6.283185307179586 \\
      --lz 3.141592653589793                       # a 2 pi x 2 x pi box
"""

import argparse
import math
import pathlib
import re
import sys

DNS_RE_TAU = 182.088     # Lee & Moser (2015), u_tau delta / nu at Re_b = 2857


def wall_coordinates(cells, gamma):
    """Vertex coordinates of the coarse mesh in y, from -1 to 1."""
    return [math.tanh(gamma * (-1 + 2 * j / cells)) / math.tanh(gamma)
            for j in range(cells + 1)]


def grid_arguments(cells, gamma, lx, lz):
    """Arguments of subdivided_hyper_rectangle with step sizes: the cell
    sizes in x; in y; in z, the two corners, colorize (boundary ids 0-1 at
    x = 0 and lx, 2-3 the walls, 4-5 at z = 0 and lz)."""
    nx, ny, nz = cells
    y = wall_coordinates(ny, gamma)
    steps = [[lx / nx] * nx,
             [y[j + 1] - y[j] for j in range(ny)],
             [lz / nz] * nz]
    lists = ";".join(",".join(f"{s:.17g}" for s in d) for d in steps)
    return f"{lists} : 0, -1, 0 : {lx:.17g}, 1, {lz:.17g} : true"


def resolution_table(cells, gamma, lx, lz, degree, refinements):
    """Lines of the wall resolution table, in DNS wall units."""
    nx, ny, nz = cells
    y = wall_coordinates(ny, gamma)
    lines = [f"Q{degree} nodes in the wall units of Re_tau = {DNS_RE_TAU:g} "
             f"(node spacing = cell size / {degree}):",
             "  refinement  cells          DoFs (u, p)  first node y+  "
             "center cell y+  dx+    dz+"]
    for r in range(refinements + 1):
        scale = 2 ** r
        n = [c * scale for c in cells]
        dofs = 4 * (degree * n[0]) * (degree * n[1] + 1) * (degree * n[2])
        first = (y[1] - y[0]) / scale / degree * DNS_RE_TAU
        center = (y[ny // 2 + 1] - y[ny // 2]) / scale * DNS_RE_TAU
        dx = lx / n[0] / degree * DNS_RE_TAU
        dz = lz / n[2] / degree * DNS_RE_TAU
        lines.append(f"  {r:>10}  {'x'.join(map(str, n)):<13}  {dofs:>11,}  "
                     f"{first:>13.2f}  {center:>14.1f}  {dx:>5.1f}  {dz:>5.1f}")
    return lines


def write_prm(path, arguments):
    text = path.read_text()
    pattern = re.compile(r"^(\s*set\s+grid arguments\s*=).*$", re.M)
    if len(pattern.findall(text)) != 1:
        sys.exit(f"error: {path} needs exactly one 'set grid arguments' line")
    path.write_text(pattern.sub(lambda m: f"{m.group(1)} {arguments}", text))


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--cells", type=int, nargs=3, default=[32, 16, 32],
                        metavar=("NX", "NY", "NZ"),
                        help="cells of the coarse mesh (default: 32 16 32)")
    parser.add_argument("--gamma", type=float, default=2.0,
                        help="stretching in y (default: 2)")
    parser.add_argument("--lx", type=float, default=4 * math.pi,
                        help="streamwise length (default: 4 pi)")
    parser.add_argument("--lz", type=float, default=2 * math.pi,
                        help="spanwise length (default: 2 pi)")
    parser.add_argument("--degree", type=int, default=2,
                        help="element degree, for the table (default: 2)")
    parser.add_argument("--refinements", type=int, default=2,
                        help="refinements in the table (default: 2)")
    parser.add_argument("--prm", type=pathlib.Path,
                        help="write the grid arguments into this .prm")
    args = parser.parse_args()
    if args.cells[1] % 2:
        sys.exit("error: NY must be even (a vertex on the center plane)")

    arguments = grid_arguments(args.cells, args.gamma, args.lx, args.lz)
    print("set grid arguments = " + arguments)
    print()
    print("\n".join(resolution_table(args.cells, args.gamma, args.lx, args.lz,
                                     args.degree, args.refinements)))
    if args.prm:
        write_prm(args.prm, arguments)
        print(f"\nwrote the grid arguments into {args.prm}")


if __name__ == "__main__":
    main()
