# SPDX-FileCopyrightText: Copyright (c) 2026
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

"""
Plot Q1 and Q2 MMS spatial convergence for the solid-phase Euler-Euler solver.

Q1 files:

    solid_mms_error_with_drag_r3.dat
    solid_mms_error_with_drag_r4.dat
    solid_mms_error_with_drag_r5.dat
    solid_mms_error_with_drag_r6.dat

Q2 files:

    solid_mms_error_com_r3_fe2.dat
    solid_mms_error_com_r4_fe2.dat
    solid_mms_error_com_r5_fe2.dat
    solid_mms_error_com_r6_fe2.dat

Each file must contain columns:

    time step alpha_L2 alpha_Linf velocity_L2 velocity_Linf

The script extracts the final row of each file, computes the mesh size h,
computes the observed order of convergence, and plots Q1 and Q2 curves
together with shaded observed-order triangles.
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

FILE_GROUPS = {
    "Q1": {
        "pattern": "solid_mms_error_with_drag_r*.dat",
        "marker": "o",
        "linestyle": "-",
    },
    "Q2": {
        "pattern": "solid_mms_error_com_r*_fe2.dat",
        "marker": "s",
        "linestyle": "-",
    },
}

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

        solid_mms_error_with_drag_r3.dat
        solid_mms_error_com_r3_fe2.dat
    """

    match = re.search(r"_r(\d+)", Path(filename).name)

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


def load_group_errors(group_name, group_settings):
    """
    Load all files for one group, for example Q1 or Q2.
    """

    script_dir = Path(__file__).resolve().parent
    files = sorted(script_dir.glob(group_settings["pattern"]))

    if len(files) == 0:
        raise RuntimeError(
            f"No files found for {group_name} with pattern "
            f"{group_settings['pattern']} in {script_dir}"
        )

    data = []

    for filename in files:
        refinement = extract_refinement_level(filename)
        row = load_error_file(filename)

        h = DOMAIN_LENGTH / (2.0 ** refinement)

        data.append({
            "group": group_name,
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


def add_order_triangle(ax, h, error, order_value):
    """
    Add a small shaded order triangle attached to the fitted convergence line.
    The sloped edge follows the fitted line, and the shaded region lies below it.
    """

    if len(h) < 2:
        return

    log_h = np.log10(h)
    log_error = np.log10(error)

    a, b = fit_log(log_h, log_error)

    # Put the triangle near the finest mesh.
    x_left = h[-1] * 1.20
    x_right = h[-1] * 1.65

    # Fitted-line values at x_left and x_right.
    y_left = 10.0 ** (a * np.log10(x_left) + b)
    y_right = 10.0 ** (a * np.log10(x_right) + b)

    x_triangle = np.array([x_left, x_right, x_right])
    y_triangle = np.array([y_left, y_left, y_right])

    ax.fill(x_triangle, y_triangle, color="gray", alpha=0.25)

    x_text = np.sqrt(x_left * x_right)
    y_text = np.sqrt(y_left * y_right) * 0.85

    ax.text(
        x_text,
        y_text,
        f"{order_value:.{n_round}f}",
        fontsize=order_font,
        ha="center",
        va="center",
    )


# -------------------------------------------------------------------------
# Printing
# -------------------------------------------------------------------------

def print_group_data(group_name, data):
    """
    Print loaded files, final-row errors, and observed orders.
    """

    refinement = np.array([row["refinement"] for row in data])
    h = np.array([row["h"] for row in data])

    alpha_L2 = np.array([row["alpha_L2"] for row in data])
    velocity_L2 = np.array([row["velocity_L2"] for row in data])

    print(f"\nLoaded MMS files: {group_name}")
    print("------------------------------------------------------------")
    for row in data:
        print(
            f"r = {row['refinement']}, "
            f"h = {row['h']:.6e}, "
            f"time = {row['time']:.6e}, "
            f"step = {row['step']}, "
            f"file = {row['filename']}"
        )

    print(f"\nFinal-row errors: {group_name}")
    print("------------------------------------------------------------")
    print("refinement        h          alpha_L2       velocity_L2")
    for row in data:
        print(
            f"{row['refinement']:10d} "
            f"{row['h']:12.5e} "
            f"{row['alpha_L2']:14.6e} "
            f"{row['velocity_L2']:14.6e}"
        )

    if len(data) >= 2:
        alpha_order = compute_pairwise_order(h, alpha_L2)
        velocity_order = compute_pairwise_order(h, velocity_L2)

        alpha_mean_order = fit_log(np.log10(h), np.log10(alpha_L2))[0]
        velocity_mean_order = fit_log(np.log10(h), np.log10(velocity_L2))[0]

        print(f"\n{group_name}: alpha L2")
        print("------------------------------------------------------------")
        print(f"Mean order from log-log fit: {alpha_mean_order:.4f}")
        print("Pairwise observed orders:")
        for i in range(len(alpha_order)):
            print(
                f"  refinement {int(refinement[i])} -> "
                f"{int(refinement[i + 1])}: {alpha_order[i]:.4f}"
            )

        print(f"\n{group_name}: velocity L2")
        print("------------------------------------------------------------")
        print(f"Mean order from log-log fit: {velocity_mean_order:.4f}")
        print("Pairwise observed orders:")
        for i in range(len(velocity_order)):
            print(
                f"  refinement {int(refinement[i])} -> "
                f"{int(refinement[i + 1])}: {velocity_order[i]:.4f}"
            )


# -------------------------------------------------------------------------
# Plotting
# -------------------------------------------------------------------------

def plot_one_quantity(all_data, quantity_key, ylabel, output_name):
    """
    Plot Q1 and Q2 results for one quantity in the same figure.
    """

    fig, ax = plt.subplots(figsize=(8, 6), facecolor="#e9e9e9")
    ax.set_facecolor("#eeeeee")

    for group_name, data in all_data.items():
        refinement = np.array([row["refinement"] for row in data])
        h = np.array([row["h"] for row in data])
        error = np.array([row[quantity_key] for row in data])

        log_h = np.log10(h)
        log_error = np.log10(error)

        coeffs = fit_log(log_h, log_error)
        mean_order = coeffs[0]

        if quantity_key == "alpha_L2":
            label = rf"{group_name.lower()}_$\alpha_s$, $p={mean_order:.{n_round}f}$"
        elif quantity_key == "velocity_L2":
            label = rf"{group_name.lower()}_$\mathbf{{u}}_s$, $p={mean_order:.{n_round}f}$"
        else:
            label = rf"{group_name.lower()}, $p={mean_order:.{n_round}f}$"

        

        ax.plot(
            h,
            error,
            marker=FILE_GROUPS[group_name]["marker"],
            linestyle=FILE_GROUPS[group_name]["linestyle"],
            linewidth=2.5,
            markersize=8,
            label=label,
        )

        add_order_triangle(ax, h, error, mean_order)

    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.set_xlabel("h")
    ax.set_ylabel(ylabel)

    ax.grid(True, which="both", linestyle=":", color="gray", alpha=0.7)

    font_properties = font_manager.FontProperties(size=11)
    ax.legend(prop=font_properties, loc="lower right")

    fig.tight_layout()

    fig.savefig(output_name + ".png", dpi=300, bbox_inches="tight")

    print(f"\nSaved {output_name}.png")

    plt.show()


# -------------------------------------------------------------------------
# Main script
# -------------------------------------------------------------------------

all_data = {}

for group_name, group_settings in FILE_GROUPS.items():
    data = load_group_errors(group_name, group_settings)
    all_data[group_name] = data
    print_group_data(group_name, data)


plot_one_quantity(
    all_data,
    "velocity_L2",
    r"$\Vert e_{\mathbf{u}_s} \Vert_2$",
    "q1_q2_order_of_convergence_solid_velocity",
)

plot_one_quantity(
    all_data,
    "alpha_L2",
    r"$\Vert e_{\alpha_s} \Vert_2$",
    "q1_q2_order_of_convergence_alpha",
)