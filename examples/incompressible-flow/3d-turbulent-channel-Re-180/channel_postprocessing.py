#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

"""Postprocessing of the turbulent channel flow at Re_tau = 180.

Compares the statistics of one or more runs of channel.prm with the DNS of
Lee & Moser (2015) at the same bulk Reynolds number, Re_b = U_b delta / nu =
2857 (Re_tau = 182).

A CASE is a folder that holds the parameter file of a run (channel.prm, or
the only .prm of the folder), the outputs of the run and its standard output
in log.txt:

    mpirun -np 192 lethe-fluid-matrix-free channel.prm | tee log.txt

Two steps:

  extract    For every snapshot written at least one time unit delta / U_b
             after the start of the averaging, the x-z average of the
             time-averaged fields (mean velocity, Reynolds stresses), folded
             onto one half-channel, in <output path>/profiles_<time>.tsv.
  (default)  Extracts the missing profiles, then compares the last snapshot
             of each case with the DNS: prints Re_tau, the error on the skin
             friction coefficient C_f and the errors on the profiles, and
             writes in --out
               channel-mean-velocity.png           U+ against y+
               channel-rms-velocities.png          u'+, v'+ and w'+ (rms)
               channel-reynolds-shear-stress.png   -u'v'+

Wall units use the friction velocity of the mean forcing. The flow control
holds the bulk velocity with a body force beta, printed in log.txt at every
time step, and in the mean beta delta = tau_w = u_tau^2 (density 1). u_tau is
the mean of beta over the averaging window of the snapshot.

Lethe's Reynolds stresses are point-wise time averages of (u - <u>)^2, with
<u> the running mean at that point. The plane average adds the x-z variance
of the time-mean field, so that the stresses are relative to the plane mean
U(y), as in the DNS.

Usage:
  python3 channel_postprocessing.py                     # the run in this folder
  python3 channel_postprocessing.py coarse fine --label "64x32x64" "128x64x128"
  python3 channel_postprocessing.py extract coarse
"""

import argparse
import base64
import math
import pathlib
import re
import sys
import xml.etree.ElementTree as ET
import zlib

import numpy as np

REFERENCE = pathlib.Path(__file__).resolve().parent / "reference"
DNS_RE_TAU = 182.088          # Lee & Moser (2015), Re_tau = 180 case
DNS_U_TAU = 6.37309e-2        # u_tau / U_b of the DNS
DNS_LABEL = "Lee & Moser (2015), DNS"
SNAPSHOT_FIELDS = ("average_velocity", "reynolds_normal_stress",
                   "reynolds_shear_stress_uv")
PROFILE_COLUMNS = ("y", "U", "uu", "vv", "ww", "uv", "dU")

# Colors of the cases (those of the other Lethe examples, in this order) and
# line styles as a second encoding. The DNS is drawn in black.
COLORS = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02"]
CASE_STYLES = list(zip(COLORS, ["--", "-", "-.", ":"]))


# --------------------------------------------------------------------------
# Parameter file
# --------------------------------------------------------------------------

def read_prm(path):
    """Entries of a deal.II parameter file, {"sub/section/key": value}, with
    the subsection names and keys in lower case and single-spaced."""
    entries, stack = {}, []
    for line in path.read_text(errors="replace").splitlines():
        line = line.split("#", 1)[0].strip()
        words = line.split(None, 1)
        if not words:
            continue
        if words[0] == "subsection" and len(words) == 2:
            stack.append(" ".join(words[1].split()).lower())
        elif words[0] == "end":
            if stack:
                stack.pop()
        elif words[0] == "set" and "=" in line:
            key, value = line[3:].split("=", 1)
            entries["/".join(stack + [" ".join(key.split()).lower()])] = \
                value.strip()
    return entries


def number(prm, key, default=math.nan):
    try:
        return float(prm[key])
    except (KeyError, ValueError):
        return default


def find_prm(path):
    """The parameter file of a case: the path itself, channel.prm in the
    folder, or the only .prm of the folder."""
    if path.is_file():
        return path
    if (path / "channel.prm").is_file():
        return path / "channel.prm"
    prms = sorted(path.glob("*.prm"))
    if len(prms) != 1:
        sys.exit(f"error: no channel.prm and {len(prms)} .prm files in {path}")
    return prms[0]


def case_parameters(prm_path):
    """What the postprocessing needs from the parameter file of a run."""
    prm = read_prm(prm_path)
    arguments = [p.strip() for p in prm["mesh/grid arguments"].split(":")]
    refinement = int(number(prm, "mesh/initial refinement", 0))
    cells = [len(d.split(",")) * 2 ** refinement
             for d in arguments[0].split(";")]
    lower = [float(c) for c in arguments[1].split(",")]
    upper = [float(c) for c in arguments[2].split(",")]
    delta = (upper[1] - lower[1]) / 2
    bulk = number(prm, "flow control/average velocity", 1.0)
    return {
        "nu": number(prm, "physical properties/fluid 0/kinematic viscosity"),
        "average_start": number(
            prm, "post-processing/initial time for average velocity", 0.0),
        "time_unit": delta / bulk,
        "bulk": bulk,
        "delta": delta,
        "cells": cells,
        "degree": int(number(prm, "fem/velocity degree", 1)),
        "output": prm_path.parent / prm.get("simulation control/output path",
                                            "./output/"),
        "output_name": prm.get("simulation control/output name", "out"),
    }


def has_statistics(params, time):
    """The snapshot written at the start of the averaging holds none."""
    return time >= params["average_start"] + params["time_unit"]


# --------------------------------------------------------------------------
# Snapshots: deal.II .vtu files, read with numpy only
# --------------------------------------------------------------------------

VTK_TYPES = {"Float32": "<f4", "Float64": "<f8", "Int8": "i1", "UInt8": "u1",
             "Int32": "<i4", "UInt32": "<u4", "Int64": "<i8", "UInt64": "<u8"}


def base64_chars(n_bytes):
    return 4 * -(-n_bytes // 3)


def decode_binary(text, dtype, compressed, header_type):
    """Values of a format="binary" DataArray: base64 of a header (block
    count, sizes) and of the zlib-compressed blocks, encoded separately."""
    raw = "".join(text.split()).encode()
    size = np.dtype(header_type).itemsize
    if not compressed:
        data = base64.b64decode(raw)
        n_bytes = int(np.frombuffer(data[:size], header_type)[0])
        return np.frombuffer(data[size:size + n_bytes], dtype)
    first = base64.b64decode(raw[:base64_chars(size)])
    n_blocks = int(np.frombuffer(first[:size], header_type)[0])
    header_bytes = (3 + n_blocks) * size
    header = np.frombuffer(
        base64.b64decode(raw[:base64_chars(header_bytes)])[:header_bytes],
        header_type)
    blob = base64.b64decode(raw[base64_chars(header_bytes):])
    blocks, position = [], 0
    for block_size in header[3:].astype(np.int64):
        blocks.append(zlib.decompress(blob[position:position + block_size]))
        position += block_size
    return np.frombuffer(b"".join(blocks), dtype)


def read_vtu(path, names):
    """Points and the named point fields of a .vtu file, all its pieces
    concatenated."""
    compressed, header_type = False, "<u4"
    points, fields, stack = [], {n: [] for n in names}, []
    for event, elem in ET.iterparse(path, events=("start", "end")):
        if event == "start":
            if elem.tag == "VTKFile":
                compressed = elem.get("compressor") is not None
                if elem.get("header_type") == "UInt64":
                    header_type = "<u8"
            stack.append(elem.tag)
            continue
        stack.pop()
        if elem.tag != "DataArray":
            continue
        parent = stack[-1] if stack else ""
        name = elem.get("Name")
        if parent == "Points" or (parent == "PointData" and name in names):
            dtype = VTK_TYPES[elem.get("type")]
            if elem.get("format", "ascii") == "binary":
                values = decode_binary(elem.text or "", dtype, compressed,
                                       header_type)
            else:
                values = np.array((elem.text or "").split(), dtype=dtype)
            values = values.astype(np.float64)
            components = int(elem.get("NumberOfComponents", "1"))
            if components > 1:
                values = values.reshape(-1, components)
            (points if parent == "Points" else fields[name]).append(values)
        elem.clear()
    missing = [n for n in names if len(fields[n]) != len(points)]
    if not points or missing:
        sys.exit(f"error: {path} has no " + ", ".join(missing or ["points"]))
    return (np.concatenate(points),
            {n: np.concatenate(v) for n, v in fields.items()})


def read_pvtu(path, names):
    """Points and fields of all the pieces of a .pvtu snapshot."""
    data = [read_vtu(path.parent / piece.get("Source"), names)
            for piece in ET.parse(path).getroot().iter("Piece")]
    return (np.concatenate([d[0] for d in data]),
            {n: np.concatenate([d[1][n] for d in data]) for n in names})


def snapshots(params):
    """(time, .pvtu path) of the snapshots of a run, from its .pvd."""
    pvd = params["output"] / f"{params['output_name']}.pvd"
    if not pvd.is_file():
        return []
    return [(float(d.get("timestep")), pvd.parent / d.get("file"))
            for d in ET.parse(pvd).getroot().iter("DataSet")
            if (pvd.parent / d.get("file")).is_file()]


def cluster(values, tolerance):
    """Labels of values grouped when closer than tolerance, and the centers."""
    order = np.argsort(values, kind="stable")
    new = np.concatenate(([True], np.diff(values[order]) > tolerance))
    labels = np.empty(len(values), dtype=np.int64)
    labels[order] = np.cumsum(new) - 1
    centers = np.bincount(labels, weights=values) / np.bincount(labels)
    return labels, centers


def plane_profiles(points, fields):
    """x-z averages of the time-averaged fields, folded onto one half: rows of
    (distance to the wall, U, u'u', v'v', w'w', u'v', U lower - U upper)."""
    length = np.ptp(points, axis=0)
    ix, x = cluster(points[:, 0], 1e-4 * length[0])
    iy, y = cluster(points[:, 1], 1e-5 * length[1])
    iz, z = cluster(points[:, 2], 1e-4 * length[2])
    # A snapshot repeats the nodes shared by cells, and the periodic images
    # at the end of the box in x and z: keep each node once.
    inside = ((x[ix] < x[-1] - 1e-3 * length[0]) &
              (z[iz] < z[-1] - 1e-3 * length[2]))
    key = (iy * (len(x) + 1) + ix) * (len(z) + 1) + iz
    _, first = np.unique(np.where(inside, key, -1), return_index=True)
    first = first[inside[first]]
    level = iy[first]
    count = np.bincount(level, minlength=len(y))

    def mean(values):
        return np.bincount(level, weights=values, minlength=len(y)) / count

    velocity = fields["average_velocity"][first]
    normal = fields["reynolds_normal_stress"][first]
    mean_u = [mean(velocity[:, i]) for i in range(3)]
    dev = [velocity[:, i] - mean_u[i][level] for i in range(3)]
    U = mean_u[0]
    uu = mean(normal[:, 0]) + mean(dev[0] ** 2)
    vv = mean(normal[:, 1]) + mean(dev[1] ** 2)
    ww = mean(normal[:, 2]) + mean(dev[2] ** 2)
    uv = mean(fields["reynolds_shear_stress_uv"][first]) + mean(dev[0] * dev[1])

    n = len(y)
    lower = np.arange((n + 1) // 2)
    upper = n - 1 - lower
    return np.column_stack([
        (y[lower] - y[0] + y[-1] - y[upper]) / 2,
        (U[lower] + U[upper]) / 2,
        (uu[lower] + uu[upper]) / 2,
        (vv[lower] + vv[upper]) / 2,
        (ww[lower] + ww[upper]) / 2,
        (uv[lower] - uv[upper]) / 2,     # u'v' changes sign across the center
        U[lower] - U[upper],
    ]), int(count.max())


def profile_path(params, time):
    return params["output"] / f"profiles_{time:g}.tsv"


def extract(params):
    """Write the profiles of the snapshots that hold statistics and have none
    yet; returns the number written."""
    written = 0
    for time, pvtu in snapshots(params):
        out = profile_path(params, time)
        if not has_statistics(params, time) or out.is_file():
            continue
        rows, nodes = plane_profiles(*read_pvtu(pvtu, SNAPSHOT_FIELDS))
        header = (f"# time {time:g}  average_start {params['average_start']:g}"
                  f"  nodes_per_plane {nodes}\n" + "\t".join(PROFILE_COLUMNS))
        np.savetxt(out, rows, delimiter="\t", header=header, comments="",
                   fmt="%.9e")
        print(f"  wrote {out}")
        written += 1
    return written


def last_profiles(params):
    """(time, profiles) of the last extracted snapshot, or None."""
    found = []
    for path in params["output"].glob("profiles_*.tsv"):
        try:
            time = float(path.stem.split("_", 1)[1])
        except ValueError:
            continue
        if has_statistics(params, time):
            found.append((time, path))
    if not found:
        return None
    time, path = max(found)
    data = np.loadtxt(path, skiprows=2, ndmin=2)
    return time, {name: data[:, i] for i, name in enumerate(PROFILE_COLUMNS)}


# --------------------------------------------------------------------------
# Statistics
# --------------------------------------------------------------------------

def read_beta(path):
    """(time, beta) of every time step, from the flow-control summary that
    the log prints with verbosity = verbose."""
    times, betas, time = [], [], None
    for line in path.read_text(errors="replace").splitlines():
        if line.startswith("Transient iteration:"):
            match = re.search(r"Time:\s*(\S+)", line)
            time = float(match.group(1)) if match else None
        elif line.startswith("Beta force:") and time is not None:
            times.append(time)
            betas.append(float(line.split(":")[1]))
            time = None
    return np.array(times), np.array(betas)


def load_dns():
    mean = np.loadtxt(REFERENCE / "LM_Channel_0180_mean_prof.dat",
                      comments="%", ndmin=2)
    fluc = np.loadtxt(REFERENCE / "LM_Channel_0180_vel_fluc_prof.dat",
                      comments="%", ndmin=2)
    return {"yp": mean[:, 1], "U": mean[:, 2], "u": np.sqrt(fluc[:, 2]),
            "v": np.sqrt(fluc[:, 3]), "w": np.sqrt(fluc[:, 4]),
            "uv": -fluc[:, 5]}


def wall_units(profiles, u_tau, nu):
    """Profiles in wall units: y+, U+, rms+, -u'v'+."""
    return {"yp": profiles["y"] * u_tau / nu,
            "U": profiles["U"] / u_tau,
            "u": np.sqrt(np.maximum(profiles["uu"], 0)) / u_tau,
            "v": np.sqrt(np.maximum(profiles["vv"], 0)) / u_tau,
            "w": np.sqrt(np.maximum(profiles["ww"], 0)) / u_tau,
            "uv": -profiles["uv"] / u_tau ** 2}


def analyze(path, label, dns):
    """Statistics of the last snapshot of a case against the DNS."""
    prm = find_prm(path)
    params = case_parameters(prm)
    extract(params)
    last = last_profiles(params)
    log = prm.parent / "log.txt"
    if last is None or not log.is_file():
        sys.exit(f"error: {path} needs log.txt and a snapshot written after "
                 f"t = {params['average_start'] + params['time_unit']:g}")
    time, profiles = last
    t, beta = read_beta(log)
    window = (t > params["average_start"]) & (t <= time + 1e-9)
    u_tau = math.sqrt(np.mean(beta[window]) * params["delta"])
    plus = wall_units(profiles, u_tau, params["nu"])

    # Relative L2 error of U+(y+), over the y+ range shared with the DNS.
    mask = plus["yp"] <= min(plus["yp"][-1], dns["yp"][-1])
    yp = plus["yp"][mask]
    reference = np.interp(yp, dns["yp"], dns["U"])
    weights = np.gradient(yp)
    u_error = math.sqrt(np.sum(weights * (plus["U"][mask] - reference) ** 2)
                        / np.sum(weights * reference ** 2))
    cells = "x".join(map(str, params["cells"]))
    return {
        "label": label or (" × ".join(map(str, params["cells"])) +
                           f" Q{params['degree']} cells"),
        "cells": cells,
        "window": (params["average_start"], time),
        "plus": plus,
        "re_tau": u_tau * params["delta"] / params["nu"],
        "cf": (u_tau / params["bulk"] / DNS_U_TAU) ** 2 - 1,
        "uc": plus["U"][-1] / dns["U"][-1] - 1,
        "u_l2": u_error,
        **{f"{k}_peak": np.max(plus[k]) / np.max(dns[k]) - 1
           for k in ("u", "v", "w", "uv")},
    }


def print_table(results):
    columns = [("Case", "label", "{}"), ("Cells", "cells", "{}"),
               ("Re_tau", "re_tau", "{:.1f}"), ("C_f", "cf", "{:+.1%}"),
               ("U_c+", "uc", "{:+.1%}"), ("U+ (L2)", "u_l2", "{:.1%}"),
               ("u'+ peak", "u_peak", "{:+.1%}"),
               ("v'+ peak", "v_peak", "{:+.1%}"),
               ("w'+ peak", "w_peak", "{:+.1%}"),
               ("-u'v'+ peak", "uv_peak", "{:+.1%}")]
    rows = [[h for h, _, _ in columns]] + [
        [form.format(r[key]) for _, key, form in columns] for r in results]
    widths = [max(len(row[i]) for row in rows) for i in range(len(columns))]
    print(f"Errors against the DNS (Re_tau = {DNS_RE_TAU}); statistics "
          "averaged from t = {:g} to {:g}".format(*results[0]["window"]))
    for row in rows:
        print("  ".join(c.rjust(w) for c, w in zip(row, widths)))


# --------------------------------------------------------------------------
# Figures
# --------------------------------------------------------------------------

def pyplot():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rc("font", weight="normal", size=13)
    plt.rcParams.update({"lines.linewidth": 2, "savefig.dpi": 300,
                         "savefig.bbox": "tight"})
    return plt


def figure(plt, results, dns, panels, path):
    """One panel per profile, against y+ on a log axis, and one legend under
    the panels."""
    fig, axes = plt.subplots(1, len(panels),
                             figsize=(max(4.6 * len(panels), 8), 4.2),
                             squeeze=False)
    for ax, (key, label) in zip(axes.flat, panels):
        # The DNS is drawn above the cases, so that it stays visible where
        # they agree with it.
        ax.plot(dns["yp"], dns[key], color="k", lw=1.5, label=DNS_LABEL,
                zorder=3)
        for result, (color, style) in zip(results, CASE_STYLES):
            plus = result["plus"]
            wall = plus["yp"] > 0     # y+ = 0 is off a log axis
            ax.plot(plus["yp"][wall], plus[key][wall], color=color, ls=style,
                    label=result["label"])
        ax.set_xscale("log")
        ax.set_xlim(0.5, 200)
        ax.set_xlabel(r"$y^+$")
        ax.set_ylabel(label)
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=len(labels),
               bbox_to_anchor=(0.5, 0.0))
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    print(f"  wrote {path}")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__.splitlines()[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="\n".join(__doc__.splitlines()[1:]))
    parser.add_argument("cases", type=pathlib.Path, nargs="*", metavar="CASE",
                        help="case folders or parameter files (default: .); "
                        "'extract' first only extracts the profiles")
    parser.add_argument("--label", nargs="+", default=[],
                        help="legend label of each case (default: its mesh)")
    parser.add_argument("--out", type=pathlib.Path, default=pathlib.Path("."),
                        help="folder of the figures (default: .)")
    args = parser.parse_args()

    if args.cases and str(args.cases[0]) == "extract":
        for case in args.cases[1:] or [pathlib.Path(".")]:
            print(f"{case}: {extract(case_parameters(find_prm(case)))} "
                  "profiles written")
        return
    cases = args.cases or [pathlib.Path(".")]
    if len(cases) > len(CASE_STYLES):
        sys.exit(f"error: at most {len(CASE_STYLES)} cases")
    if args.label and len(args.label) != len(cases):
        sys.exit("error: give one --label per case")
    labels = args.label or [None] * len(cases)

    dns = load_dns()
    results = [analyze(c, l, dns) for c, l in zip(cases, labels)]
    print_table(results)

    plt = pyplot()
    args.out.mkdir(parents=True, exist_ok=True)
    figure(plt, results, dns, [("U", r"$U^+$")],
           args.out / "channel-mean-velocity.png")
    figure(plt, results, dns, [("u", r"$u'^+_\mathrm{rms}$"),
                               ("v", r"$v'^+_\mathrm{rms}$"),
                               ("w", r"$w'^+_\mathrm{rms}$")],
           args.out / "channel-rms-velocities.png")
    figure(plt, results, dns, [("uv", r"$-\overline{u'v'}^+$")],
           args.out / "channel-reynolds-shear-stress.png")


if __name__ == "__main__":
    main()
