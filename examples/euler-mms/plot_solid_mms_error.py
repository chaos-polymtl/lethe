# SPDX-FileCopyrightText: Copyright (c) 2026
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

"""
Plot MMS spatial convergence for the solid-phase Euler-Euler solver.

This script reads files named

    solid_mms_error_r3.dat
    solid_mms_error_r4.dat
    solid_mms_error_r5.dat
    solid_mms_error_r6.dat

Each file must contain columns:

    time step alpha_L2 alpha_Linf velocity_L2 velocity_Linf

The script extracts the final row of each file, computes the mesh size h,
computes the observed order of convergence, and plots Lethe-style log-log
error curves with a shaded order triangle.
"""

import re
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import font_manager


# -------------------------------------------------------------------------
# User settings
# -------------------------------------------------------------------------

DOMAIN_LENGTH = 2.0  
FILE_PATTERN = "solid_mms_error_r*.dat"

markers = ["o", "s", "^", "D", "x", "*"]
order_font = 10
n_round = 2


# -------------------------------------------------------------------------
# Plot style
# -------------------------------------------------------------------------

plt.rcParams.update({
    "figure.facecolor": "#e9e9e9",
    "axes.facecolor": "#eeeeee",
    "axes.labelsize": 20,
    "axes.titlesize": 22,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 12,
    "text.usetex": False,
})


# -------------------------------------------------------------------------
# Data reading
# -------------------------------------------------------------------------

def extract_refinement_level(filename):
    """
    Extract refinement level from names such as:

        solid_mms_error_r3.dat
        solid_mms_error_r4.dat
    """

    match = re.search(r"_r(\d+)\.dat$", str(filename))

    if match is None:
        raise ValueError(f"Could not extract refinement level from {filename}")

    return int(match.group(1))


def load_error_file(filename):
    """
    Read one MMS error file and return the final-row errors.

    Columns:
        time step alpha_L2 alpha_Linf velocity_L2 velocity_Linf
    """

    data = np.loadtxt(filename, skiprows=1)

    if data.ndim == 1:
        last_row = data
    else:
        last_row = data[-1, :]

    return {
        "time": last_row[0],
        "step": int(last_row[1]),
        "alpha_L2": last_row[2],
        "alpha_Linf": last_row[3],
        "velocity_L2": last_row[4],
        "velocity_Linf": last_row[5],
    }


def load_all_errors():
    """
    Load all solid_mms_error_r*.dat files from the folder where this script is saved.
    """

    script_dir = Path(__file__).resolve().parent
    files = sorted(script_dir.glob(FILE_PATTERN))

    if len(files) == 0:
        raise RuntimeError(
            f"No files found with pattern {FILE_PATTERN} in {script_dir}"
        )

    data = []

    for filename in files:
        refinement = extract_refinement_level(filename)
        row = load_error_file(filename)

        h = DOMAIN_LENGTH / (2.0 ** refinement)

        data.append({
            "filename": filename.name,
            "refinement": refinement,
            "h": h,
            **row,
        })

    return sorted(data, key=lambda row: row["refinement"])


# -------------------------------------------------------------------------
# Error/order calculations
# -------------------------------------------------------------------------

def fit_log(log_h, log_error):
    """
    Fit log(error) = a log(h) + b.

    The slope a is the observed order of accuracy.
    """

    return np.polyfit(log_h, log_error, 1)


def compute_pairwise_order(h, error):
    """
    Compute observed order between consecutive refinements.
    """

    orders = []

    for i in range(len(error) - 1):
        p = np.log(error[i] / error[i + 1]) / np.log(h[i] / h[i + 1])
        orders.append(p)

    return np.array(orders)




def add_order_triangle(ax, h, error, order_value, x0=None, x1=None):
    """
    Add a shaded triangle showing the fitted order of convergence.
    """

    if x0 is None:
        x0 = h[-2]

    if x1 is None:
        x1 = h[-1]

    log_h = np.log10(h)
    log_error = np.log10(error)

    a, b = fit_log(log_h, log_error)

    y0 = 10.0 ** (a * np.log10(x0) + b)
    y1 = 10.0 ** (a * np.log10(x1) + b)

    x_triangle = np.array([x0, x1, x1])
    y_triangle = np.array([y0, y0, y1])

    ax.fill(x_triangle, y_triangle, color="gray", alpha=0.25)

    x_text = np.sqrt(x0 * x1)
    y_text = np.sqrt(y0 * y1)

    ax.text(
        x_text,
        y_text,
        f"{order_value:.{n_round}f}",
        fontsize=order_font,
        ha="center",
        va="center",
    )


# -------------------------------------------------------------------------
# Plotting
# -------------------------------------------------------------------------

def plot_one_quantity(refinement, h, error, ylabel, output_name, label, marker):
    

    log_h = np.log10(h)
    log_error = np.log10(error)

    coeffs = fit_log(log_h, log_error)
    mean_order = coeffs[0]

    pairwise_order = compute_pairwise_order(h, error)

    fig, ax = plt.subplots(figsize=(8, 6), facecolor="#e9e9e9")
    ax.set_facecolor("#eeeeee")

    ax.plot(
        h,
        error,
        marker=marker,
        linewidth=2.5,
        markersize=8,
        label=label,
    )

    add_order_triangle(
        ax,
        h,
        error,
        mean_order,
        x0=h[-2],
        x1=h[-1],
    )

    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.set_xlabel("h")
    ax.set_ylabel(ylabel)
    

    ax.grid(True, which="both", linestyle=":", color="gray", alpha=0.7)

    font_properties = font_manager.FontProperties(size=11)
    ax.legend(prop=font_properties, loc="lower right")

    fig.tight_layout()

    fig.savefig(output_name + ".png", dpi=300, bbox_inches="tight")
    

    print(f"\n{label}")
    print("------------------------------------------------------------")
    print(f"Mean order from log-log fit: {mean_order:.4f}")
    print("Pairwise observed orders:")

    for i in range(len(pairwise_order)):
        print(
            f"  refinement {int(refinement[i])} -> "
            f"{int(refinement[i + 1])}: {pairwise_order[i]:.4f}"
        )

    print(f"Saved {output_name}.png")
    
    plt.show()


# -------------------------------------------------------------------------
# Main script
# -------------------------------------------------------------------------

data = load_all_errors()

refinement = np.array([row["refinement"] for row in data])
h = np.array([row["h"] for row in data])

alpha_L2 = np.array([row["alpha_L2"] for row in data])
alpha_Linf = np.array([row["alpha_Linf"] for row in data])

velocity_L2 = np.array([row["velocity_L2"] for row in data])
velocity_Linf = np.array([row["velocity_Linf"] for row in data])


print("\nLoaded MMS files")
print("------------------------------------------------------------")
for row in data:
    print(
        f"r = {row['refinement']}, "
        f"h = {row['h']:.6e}, "
        f"time = {row['time']:.6e}, "
        f"step = {row['step']}, "
        f"file = {row['filename']}"
    )


print("\nFinal-row errors")
print("------------------------------------------------------------")
print("refinement        h          alpha_L2       velocity_L2")
for row in data:
    print(
        f"{row['refinement']:10d} "
        f"{row['h']:12.5e} "
        f"{row['alpha_L2']:14.6e} "
        f"{row['velocity_L2']:14.6e}"
    )


plot_one_quantity(
    refinement,
    h,
    velocity_L2,
    r"$\Vert e_{\mathbf{u}_s} \Vert_2$",
    "order_of_convergence_solid_velocity",
    r"$\mathbf{u}_s$, $L_2$",
    markers[0],
)

plot_one_quantity(
    refinement,
    h,
    alpha_L2,
    r"$\Vert e_{\alpha_s} \Vert_2$",
    "order_of_convergence_alpha",
    r"$\alpha_s$, $L_2$",
    markers[1],
)