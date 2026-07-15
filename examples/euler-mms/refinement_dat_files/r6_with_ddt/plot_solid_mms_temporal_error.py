

"""
Temporal MMS difference-error plot for fixed mesh size.

This script is designed for the article-style temporal verification procedure.

At fixed mesh size, the total error is

    E(dt) = C_h + g_t dt^q

where C_h is the constant spatial-error contribution.

Therefore, the raw error E(dt) may appear almost flat.
To remove the constant spatial contribution, this script plots

    |E(dt) - E(dt/r_t)|  versus  dt

For a constant temporal refinement factor r_t,

    |E(dt) - E(dt/r_t)| ~ dt^q

so the slope of the log-log plot gives the temporal order q.

Expected filenames:

    solid_mms_error_dt1e-4.dat
    solid_mms_error_dt5e-5.dat
    solid_mms_error_dt2p5e-5.dat
    solid_mms_error_dt1p25e-5.dat

Each file must contain columns:

    time step alpha_L2 alpha_Linf velocity_L2 velocity_Linf
"""

import re
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.ticker import NullFormatter, FixedLocator, FixedFormatter




# -------------------------------------------------------------------------
# User settings
# -------------------------------------------------------------------------

FILE_PATTERN = "solid_mms_error_dt*.dat"

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
    "axes.titlesize": 20,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 12,
    "text.usetex": False,
})


# -------------------------------------------------------------------------
# Filename parsing
# -------------------------------------------------------------------------

def token_to_float(token):
    """
    Convert strings such as:
        1e-4
        2p5e-5
        1p25e-5
    to float.
    """

    return float(token.replace("p", ".").replace("P", "."))


def extract_time_step(filename):
    """
    Extract dt from filenames such as:

        solid_mms_error_dt1e-4.dat
        solid_mms_error_dt5e-5.dat
        solid_mms_error_dt2p5e-5.dat
        solid_mms_error_dt2.5e-5.dat
    """

    name = Path(filename).name

    match = re.search(
        r"dt([0-9]+(?:[pP\.][0-9]+)?(?:[eE][+-]?\d+)?)",
        name,
    )

    if match is None:
        raise ValueError(f"Could not extract dt from filename: {filename}")

    return token_to_float(match.group(1))


# -------------------------------------------------------------------------
# Data reading
# -------------------------------------------------------------------------

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
    Load all temporal MMS error files from the folder where this script is saved.
    """

    script_dir = Path(__file__).resolve().parent
    files = sorted(script_dir.glob(FILE_PATTERN))

    if len(files) == 0:
        raise RuntimeError(
            f"No files found with pattern {FILE_PATTERN} in {script_dir}"
        )

    data = []

    for filename in files:
        dt = extract_time_step(filename)
        row = load_error_file(filename)

        data.append({
            "filename": filename.name,
            "dt": dt,
            **row,
        })

    # Sort from coarse to fine:
    # 1e-4, 5e-5, 2.5e-5, ...
    return sorted(data, key=lambda row: row["dt"], reverse=True)


# -------------------------------------------------------------------------
# Difference-error calculations
# -------------------------------------------------------------------------

def compute_difference_error(dt, error):
    """
    Compute temporal difference errors:

        D_i = |E(dt_i) - E(dt_{i+1})|

    The x-value is taken as the coarser time step dt_i.
    """

    dt_diff = []
    diff_error = []
    refinement_factors = []

    for i in range(len(error) - 1):
        dt_i = dt[i]
        dt_next = dt[i + 1]

        d_error = abs(error[i] - error[i + 1])

        dt_diff.append(dt_i)
        diff_error.append(d_error)
        refinement_factors.append(dt_i / dt_next)

    return (
        np.array(dt_diff),
        np.array(diff_error),
        np.array(refinement_factors),
    )


def fit_log_log_order(x, y):
    """
    Fit log(y) = q log(x) + c.
    """

    coeffs = np.polyfit(np.log10(x), np.log10(y), 1)
    return coeffs[0]


def compute_pairwise_order(x, y):
    """
    Pairwise order for the difference-error data.
    """

    orders = []

    for i in range(len(y) - 1):
        q = np.log(y[i] / y[i + 1]) / np.log(x[i] / x[i + 1])
        orders.append(q)

    return np.array(orders)


def add_order_triangle(ax, x, y, order_value):
    """
    Add a shaded order triangle attached to the fitted line.
    """

    if len(x) < 2:
        return

    a, b = np.polyfit(np.log10(x), np.log10(y), 1)

    x_left = x[-1] * 1.15
    x_right = x[-1] * 1.65

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



def format_dt_tick(value):
    """
    Format time-step values clearly for the x-axis.
    Example:
        1e-4, 5e-5, 2.5e-5, 1.25e-5
    """

    text = f"{value:.3g}"

    # Convert 0.0001 to 1e-04 style if needed
    if "e" not in text:
        text = f"{value:.2e}"

    text = text.replace("e-0", "e-")
    text = text.replace("e+0", "e+")

    return text


# -------------------------------------------------------------------------
# Plotting
# -------------------------------------------------------------------------

def plot_difference_quantity(dt, error, ylabel, output_name, label, marker):
    """
    Plot |E(dt) - E(dt/r_t)| versus dt.
    """

    dt_diff, diff_error, refinement_factors = compute_difference_error(dt, error)

    if len(diff_error) < 2:
        print(f"\n{label}")
        print("------------------------------------------------------------")
        print("Not enough points for a difference-error convergence plot.")
        print("Need at least three dt files.")
        return

    fit_order = fit_log_log_order(dt_diff, diff_error)
    pairwise_orders = compute_pairwise_order(dt_diff, diff_error)

    fig, ax = plt.subplots(figsize=(8, 6), facecolor="#e9e9e9")
    ax.set_facecolor("#eeeeee")

    ax.plot(
        dt_diff,
        diff_error,
        marker=marker,
        linewidth=2.5,
        markersize=8,
        label=label,
    )

    add_order_triangle(ax, dt_diff, diff_error, fit_order)

    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.xaxis.set_major_locator(FixedLocator(dt_diff))
    ax.xaxis.set_major_formatter(
    FixedFormatter([format_dt_tick(value) for value in dt_diff])
    )

    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.tick_params(axis="x", which="major", labelsize=11)
    ax.tick_params(axis="x", which="minor", labelbottom=False)

    ax.set_xlabel(r"$\Delta t$")
    ax.set_ylabel(ylabel)

    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.grid(True, which="both", linestyle=":", color="gray", alpha=0.7)

    font_properties = font_manager.FontProperties(size=11)
    ax.legend(prop=font_properties, loc="best")

    fig.subplots_adjust(left=0.18, right=0.96, bottom=0.15, top=0.95)

    script_dir = Path(__file__).resolve().parent
    png_file = script_dir / f"{output_name}.png"
    

    fig.savefig(png_file, dpi=300)
 

    print(f"\n{label}")
    print("------------------------------------------------------------")
    print("Difference errors:")
    for i in range(len(diff_error)):
        print(
            f"  |E({dt[i]:.6e}) - E({dt[i + 1]:.6e})| "
            f"= {diff_error[i]:.6e}, "
            f"r_t = {refinement_factors[i]:.6g}"
        )

    print(f"Order from log-log fit of difference errors: {fit_order:.6f}")

    print("Pairwise orders from difference errors:")
    for i in range(len(pairwise_orders)):
        print(
            f"  {dt_diff[i]:.6e} -> {dt_diff[i + 1]:.6e}: "
            f"{pairwise_orders[i]:.6f}"
        )

    print(f"Saved {png_file}")
    

    plt.show()


# -------------------------------------------------------------------------
# Main script
# -------------------------------------------------------------------------

data = load_all_errors()

dt = np.array([row["dt"] for row in data])

alpha_L2 = np.array([row["alpha_L2"] for row in data])
alpha_Linf = np.array([row["alpha_Linf"] for row in data])

velocity_L2 = np.array([row["velocity_L2"] for row in data])
velocity_Linf = np.array([row["velocity_Linf"] for row in data])


print("\nLoaded temporal MMS files")
print("------------------------------------------------------------")
for row in data:
    print(
        f"dt = {row['dt']:.6e}, "
        f"time = {row['time']:.6e}, "
        f"step = {row['step']}, "
        f"file = {row['filename']}"
    )


print("\nFinal-row total errors")
print("------------------------------------------------------------")
print("      dt         alpha_L2       alpha_Linf      velocity_L2    velocity_Linf")
for row in data:
    print(
        f"{row['dt']:12.5e} "
        f"{row['alpha_L2']:14.6e} "
        f"{row['alpha_Linf']:14.6e} "
        f"{row['velocity_L2']:14.6e} "
        f"{row['velocity_Linf']:14.6e}"
    )


plot_difference_quantity(
    dt,
    velocity_L2,
    r"$\left|E_{\mathbf{u}_s}(\Delta t)-E_{\mathbf{u}_s}(\Delta t/2)\right|$",
    "temporal_difference_solid_velocity_L2",
    r"$\mathbf{u}_s$, temporal difference",
    markers[0],
)

plot_difference_quantity(
    dt,
    alpha_L2,
    r"$\left|E_{\alpha_s}(\Delta t)-E_{\alpha_s}(\Delta t/2)\right|$",
    "temporal_difference_alpha_L2",
    r"$\alpha_s$, temporal difference",
    markers[1],
)