

import re
import glob
import math
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


DOMAIN_LENGTH = 2.0   # For Omega = [-1,1] x [-1,1]
EXPECTED_ORDER = 2.0  # Case C with Q1 + BDF1 + dt ~ h^2


def parse_float_from_filename_dt(dt_string):
    """
    Converts filename dt strings such as:
      4e-3
      1e-3
      2.5e-4
      6.25e-5
    into floats.
    """
    return float(dt_string)


def parse_filename(filename):
    """
    Expected filename format:
      solid_mms_error_com_r3_dt4e-3.dat
      solid_mms_error_com_r4_dt1e-3.dat
      solid_mms_error_com_r5_dt2.5e-4.dat
      solid_mms_error_com_r6_dt6.25e-5.dat
    """

    name = Path(filename).stem

    match = re.search(r"_r(\d+)_dt([0-9.+\-eE]+)$", name)

    if not match:
        raise ValueError(f"Could not parse refinement and dt from filename: {filename}")

    refinement = int(match.group(1))
    dt = parse_float_from_filename_dt(match.group(2))

    h = DOMAIN_LENGTH / (2 ** refinement)

    return refinement, h, dt


def parse_error_file(filename):
    """
    Reads alpha L2 and velocity L2 errors from either:

    1. log-style files:
         alpha L2       = ...
         velocity L2    = ...

    2. table-style .dat files, for example:
         time alpha_L2 alpha_Linf velocity_L2 velocity_Linf
         0.1  1e-5     2e-5      3e-5        4e-5

    The last valid data row is used.
    """

    path = Path(filename)
    text = path.read_text(errors="ignore")

    # ------------------------------------------------------------
    # First try log-style parsing
    # ------------------------------------------------------------
    alpha_matches = re.findall(
        r"alpha\s+L2\s*=\s*([+-]?(?:\d+\.\d*|\d*\.\d+|\d+)(?:[eE][+-]?\d+)?)",
        text,
    )

    velocity_matches = re.findall(
        r"velocity\s+L2\s*=\s*([+-]?(?:\d+\.\d*|\d*\.\d+|\d+)(?:[eE][+-]?\d+)?)",
        text,
    )

    if alpha_matches and velocity_matches:
        alpha_l2 = float(alpha_matches[-1])
        velocity_l2 = float(velocity_matches[-1])
        return alpha_l2, velocity_l2

    # ------------------------------------------------------------
    # Then try table-style parsing
    # ------------------------------------------------------------
    lines = text.splitlines()

    header = None
    numeric_rows = []

    number_pattern = re.compile(
        r"^[\s,]*[+-]?(?:\d+\.\d*|\d*\.\d+|\d+)(?:[eE][+-]?\d+)?"
    )

    for line in lines:
        clean = line.strip()

        if not clean:
            continue

        # Remove common comment symbols but keep possible header text
        possible_header = clean.lstrip("#").strip()

        # Check if line starts with a number
        if number_pattern.match(clean):
            parts = re.split(r"[\s,]+", clean)

            try:
                values = [float(p) for p in parts if p != ""]
            except ValueError:
                continue

            if len(values) >= 2:
                numeric_rows.append(values)
        else:
            # Store last non-numeric line as possible header
            header = possible_header

    if not numeric_rows:
        raise ValueError(f"No numeric data rows found in file: {filename}")

    last_row = numeric_rows[-1]

    # ------------------------------------------------------------
    # If header exists, try to identify columns by name
    # ------------------------------------------------------------
    if header is not None:
        header_parts = re.split(r"[\s,]+", header.strip())

        normalized = [
            h.lower()
            .replace("-", "_")
            .replace("(", "")
            .replace(")", "")
            .replace("[", "")
            .replace("]", "")
            for h in header_parts
        ]

        alpha_index = None
        velocity_index = None

        for i, name in enumerate(normalized):
            if name in ["alpha_l2", "alpha_l2_error", "l2_alpha", "l2_alpha_s"]:
                alpha_index = i

            if name in [
                "velocity_l2",
                "velocity_l2_error",
                "u_l2",
                "us_l2",
                "u_s_l2",
                "l2_velocity",
                "l2_u",
                "l2_us",
            ]:
                velocity_index = i

        if alpha_index is not None and velocity_index is not None:
            if alpha_index < len(last_row) and velocity_index < len(last_row):
                return last_row[alpha_index], last_row[velocity_index]

    # ------------------------------------------------------------
    # Fallback assumptions
    # ------------------------------------------------------------
    # Common MMS table format:
    # time alpha_L2 alpha_Linf velocity_L2 velocity_Linf
    if len(last_row) >= 5:
        alpha_l2 = last_row[1]
        velocity_l2 = last_row[3]
        return alpha_l2, velocity_l2

    # Possible format:
    # alpha_L2 alpha_Linf velocity_L2 velocity_Linf
    if len(last_row) == 4:
        alpha_l2 = last_row[0]
        velocity_l2 = last_row[2]
        return alpha_l2, velocity_l2

    # Possible format:
    # alpha_L2 velocity_L2
    if len(last_row) == 2:
        alpha_l2 = last_row[0]
        velocity_l2 = last_row[1]
        return alpha_l2, velocity_l2

    raise ValueError(
        f"Could not identify alpha L2 and velocity L2 columns in file: {filename}\n"
        f"Last numeric row found was: {last_row}\n"
        f"Header found was: {header}"
    )


def compute_orders(errors, hs):
    """
    Computes observed order between consecutive mesh refinements:

      p = log(E_i / E_{i+1}) / log(h_i / h_{i+1})
    """

    orders = [np.nan]

    for i in range(1, len(errors)):
        e_coarse = errors[i - 1]
        e_fine = errors[i]
        h_coarse = hs[i - 1]
        h_fine = hs[i]

        p = math.log(e_coarse / e_fine) / math.log(h_coarse / h_fine)
        orders.append(p)

    return orders


def make_reference_line(hs, errors, order):
    """
    Creates an O(h^order) reference line anchored at the first error value.
    """

    h0 = hs[0]
    e0 = errors[0]

    return np.array([e0 * (h / h0) ** order for h in hs])


def main():
    files = sorted(glob.glob("solid_mms_error_com_r*_dt*.dat"))

    if not files:
        raise FileNotFoundError(
            "No files found matching: solid_mms_error_com_r*_dt*.dat"
        )

    data = []

    for filename in files:
        refinement, h, dt = parse_filename(filename)
        alpha_l2, velocity_l2 = parse_error_file(filename)

        data.append(
            {
                "file": filename,
                "r": refinement,
                "h": h,
                "dt": dt,
                "alpha_l2": alpha_l2,
                "velocity_l2": velocity_l2,
            }
        )

    data.sort(key=lambda row: row["r"])

    refinements = np.array([row["r"] for row in data])
    hs = np.array([row["h"] for row in data])
    dts = np.array([row["dt"] for row in data])
    alpha_errors = np.array([row["alpha_l2"] for row in data])
    velocity_errors = np.array([row["velocity_l2"] for row in data])

    alpha_orders = compute_orders(alpha_errors, hs)
    velocity_orders = compute_orders(velocity_errors, hs)

    print("\nCase C MMS convergence results")
    print("-" * 95)
    print(
        f"{'r':>4} {'h':>14} {'dt':>14} "
        f"{'alpha L2':>16} {'alpha order':>14} "
        f"{'velocity L2':>16} {'velocity order':>16}"
    )
    print("-" * 95)

    for i, row in enumerate(data):
        alpha_order_text = "--" if i == 0 else f"{alpha_orders[i]:.4f}"
        velocity_order_text = "--" if i == 0 else f"{velocity_orders[i]:.4f}"

        print(
            f"{row['r']:4d} "
            f"{row['h']:14.6e} "
            f"{row['dt']:14.6e} "
            f"{row['alpha_l2']:16.6e} "
            f"{alpha_order_text:>14} "
            f"{row['velocity_l2']:16.6e} "
            f"{velocity_order_text:>16}"
        )

    print("-" * 95)

    # Reference lines
    alpha_ref = make_reference_line(hs, alpha_errors, EXPECTED_ORDER)
    velocity_ref = make_reference_line(hs, velocity_errors, EXPECTED_ORDER)

    # Plot
    plt.figure(figsize=(7.0, 5.2))

    plt.loglog(
        hs,
        alpha_errors,
        marker="o",
        linewidth=2,
        label=r"$\alpha_s$ $L_2$ error",
    )

    plt.loglog(
        hs,
        velocity_errors,
        marker="s",
        linewidth=2,
        label=r"$\mathbf{u}_s$ $L_2$ error",
    )

    plt.loglog(
        hs,
        alpha_ref,
        linestyle="--",
        linewidth=1.5,
        label=r"$O(h^2)$ reference for $\alpha_s$",
    )

    plt.loglog(
        hs,
        velocity_ref,
        linestyle=":",
        linewidth=1.8,
        label=r"$O(h^2)$ reference for $\mathbf{u}_s$",
    )

    plt.gca().invert_xaxis()

    plt.xlabel(r"Mesh size, $h$")
    plt.ylabel(r"$L_2$ error")
    plt.title(r"Case C: fully discrete MMS convergence")
    plt.grid(True, which="both", linestyle="--", linewidth=0.5)
    plt.legend()
    plt.tight_layout()

    output_name = "case_C_mms_convergence.png"
    plt.savefig(output_name, dpi=300)

    print(f"\nSaved plot: {output_name}\n")


if __name__ == "__main__":
    main()