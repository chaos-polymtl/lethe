# SPDX-FileCopyrightText: Copyright (c) 2025-2026 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#############################################################################
# Import Libraries
#############################################################################
import argparse
import numpy as np

import matplotlib.pyplot as plt

#############################################################################
# Setup plot parameters
#############################################################################

# Marker by bed type
marker_map = {
    "MB PCM": "o",
    "MB QCM": "s",
    "MF QCM": "D"
}

# Color by drag/coupling scheme
color_map = {
    "Explicit": "#1b9e77",   
    "Semi-Implicit": "#d95f02",  
    "Implicit": "#7570b3"   
}


#############################################################################
# Define argument parser
#############################################################################
parser = argparse.ArgumentParser(
    description='Arguments for the post-processing of the cylindrical fluidized bed example'
)
parser.add_argument(
    "--validate", action="store_true",
    help="Launches the script in validation mode. This will log the content of the graph and prevent the display of figures",
    default=False)
parser.add_argument(
    "--solver", type=str, nargs="+", choices=["mb", "mf"],
    help="Solver type(s): mb or mf"
)
parser.add_argument(
    "--filter", type=str, nargs="+", choices=["pcm", "qcm"],
    help="Filter type(s): pcm or qcm"
)
parser.add_argument(
    "--drag", type=str, nargs="+",
    choices=["explicit", "semi-implicit", "implicit"],
    help="Drag coupling scheme(s)"
)
parser.add_argument( "--output-suffix", "-s", type=str, help="Suffix for the plot name")

args, leftovers = parser.parse_known_args()

solver_given = args.solver is not None
filter_given = args.filter is not None
drag_given   = args.drag   is not None

if args.output_suffix and not (solver_given and filter_given and drag_given):
    raise ValueError(
        "If you use --output-suffix, you must also provide --solver, --filter, and --drag."
    )
if (solver_given or filter_given or drag_given) and not (solver_given and filter_given and drag_given):
    raise ValueError("You must provide --solver, --filter, and --drag together.")
    

#############################################################################
# Helper functions
#############################################################################
# === Binning function ===
def bin_and_average(time, pressure, bins, averaging_window):
    indices = np.digitize(time, bins) - 1
    mean_p, std_p = [], []
    for i in range(len(bins)-1):
        mask = (indices == i) & (time > (bins[i+1] - averaging_window))
        if np.any(mask):
            mean_p.append(pressure[mask].mean())
            std_p.append(pressure[mask].std())
        else:
            mean_p.append(np.nan)
            std_p.append(np.nan)
    return np.array(mean_p), np.array(std_p)

# === Plotting function ===
def plot_pressure(datasets, Re, dataset_keys, subscript, validate):
    for i, key in enumerate(dataset_keys):
        # Determine which marker (corresponds to combination of solver and filter used) and what color (corresponds to drag scheme) to use for the plot
        parts = key.split()
        solver = " ".join(parts[:2])
        drag_coupling = parts[-1]
        
        marker = marker_map[solver]
        color = color_map[drag_coupling]
        
        plt.errorbar(Re, datasets[key]["mean"], yerr=datasets[key]["std"], fmt=marker+'-',  color=color, markerfacecolor='none', markersize=8, capsize=3, label=key)
        
    plt.plot([0, max(Re)], [delta_p_analytical]*2, 'k--', label="Fluidization Δp")
    plt.plot([Re_mf_WY_inlet, Re_mf_WY_inlet], [0, max([datasets[k]["mean"].max() for k in dataset_keys])*1.1], 'k:')
    plt.plot([Re_mf_N_inlet, Re_mf_N_inlet], [0, max([datasets[k]["mean"].max() for k in dataset_keys])*1.1], 'k-.')
      # Add annotations below the lines with arrows pointing up  
    # Fixed y for annotations
    y_annot = 80
    
    plt.annotate(
        "Wen-Yu",
        xy=(Re_mf_WY_inlet, y_annot),        
        xytext=(Re_mf_WY_inlet - 0.1*max(Re), y_annot),  
        arrowprops=dict(facecolor='black', arrowstyle="->",
                        connectionstyle="arc3,rad=0.3"),    
        horizontalalignment='right',
        verticalalignment='top',
        fontsize=12
    )
    # Noda annotation (arrow from right, curved)
    plt.annotate(
        "Noda et al.",
        xy=(Re_mf_N_inlet, y_annot),
        xytext=(Re_mf_N_inlet + 0.09*max(Re), y_annot),  
        arrowprops=dict(facecolor='black', arrowstyle="->",
                        connectionstyle="arc3,rad=-0.3"),  
        horizontalalignment='left',
        verticalalignment='top',
        fontsize=12
    )
    plt.grid(True, which='major', linestyle='--', linewidth=0.5, alpha=0.6)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlabel("Re", fontsize=14)
    plt.ylabel("$\\Delta p$ (Pa)", fontsize=14)
    plt.legend()
    if validate:
        plt.savefig(f"pressure-drop-Re-{subscript}.pdf")
        solution = np.hstack([datasets[k]["mean"].reshape(-1,1) for k in dataset_keys])
        return solution
    else:
        plt.savefig(f"pressure-drop-{subscript}.png", dpi=600, bbox_inches='tight')
        plt.show()

# Use input data to build the keys of the dataset to be used in the plot
def build_key(solver, flt, drag):
    # Pretty label used in plots
    solver_label = solver.upper()
    filter_label = flt.upper()
    drag_label = drag.title() if drag != "semi-implicit" else "Semi-Implicit"
    return f"{solver_label} {filter_label} {drag_label}"

def build_path(solver, flt, drag):
    # Matches your directory naming convention
    solver_part = "mb" if solver == "mb" else "mf"
    filter_part = flt
    drag_part = drag.replace("-", "_")

    if solver == "mb":
        return f"output_{solver_part}_{filter_part}_{drag_part}/pressure_drop.dat"
    else:
        return f"output_{solver_part}_{drag_part}/pressure_drop.dat"

#############################################################################
# Load data
#############################################################################
datasets = {}
data_files = {}

if args.solver and args.filter and args.drag:
    for solver in args.solver:
        for flt in args.filter:
            for drag in args.drag:

                # MF only supports QCM
                if solver == "mf" and flt == "pcm":
                    continue

                key = build_key(solver, flt, drag)
                path = build_path(solver, flt, drag)

                data_files[key] = path
else:
    data_files = {
        "MB PCM Explicit": "output_mb_pcm_explicit/pressure_drop.dat",
        "MB PCM Semi-Implicit": "output_mb_pcm_semi_implicit/pressure_drop.dat",
        "MB PCM Implicit": "output_mb_pcm_implicit/pressure_drop.dat",
        "MB QCM Explicit": "output_mb_qcm_explicit/pressure_drop.dat",
        "MB QCM Semi-Implicit": "output_mb_qcm_semi_implicit/pressure_drop.dat",
        "MB QCM Implicit": "output_mb_qcm_implicit/pressure_drop.dat",
        "MF QCM Explicit": "output_mf_explicit/pressure_drop.dat",
        "MF QCM Semi-Implicit": "output_mf_semi_implicit/pressure_drop.dat",
        "MF QCM Implicit": "output_mf_implicit/pressure_drop.dat"
    }

for key, path in data_files.items():
    t, p, _ = np.loadtxt(path, unpack=True, skiprows=1)
    datasets[key] = {"t": t, "p": p}

print(datasets)
#############################################################################
# Physical properties
#############################################################################    
# Column diameter
D= 0.02
# Bed cross sectional area
Ab = np.pi*(D/2)**2

# Particle properties
dp = 5e-4
Vp = 4/3*np.pi*(dp/2)**3
rho_p = 1000
Np = 200000

# Fluid properties
nu = 1e-5
rho_f = 1
mu = nu*rho_f

g = 9.81
Ar =  g * rho_f * (rho_p-rho_f)*dp**3 / mu**2
print("Ar: ", Ar)

#############################################################################
# Existing correlations for minimum fluidization velocity
#############################################################################
# Minimum fluidization velocity
Re_mf_WY = (33.7**2+0.0408*Ar)**0.5 - 33.7
U_mf_WY = Re_mf_WY * nu / dp 
Re_mf_WY_inlet = U_mf_WY * D/nu

Re_mf_N = (19.29**2+0.0276*Ar)**0.5 - 19.29
U_mf_N = Re_mf_N * nu / dp
Re_mf_N_inlet= U_mf_N * D/nu

print("Wen-Yu: Re_p =", Re_mf_WY, "U_mf =", U_mf_WY, "Re_inlet =", Re_mf_WY_inlet)
print("Noda: Re_p =", Re_mf_N, "U_mf =", U_mf_N, "Re_inlet =", Re_mf_N_inlet)

#############################################################################
# Pressure drop after fluidization
#############################################################################
delta_p_analytical = Np*Vp*(rho_p-rho_f)*g / Ab
print("Analytical pressure drop: ", delta_p_analytical)

#############################################################################
# Average pressure drop for each velocity step
#############################################################################
# Time interval for each velocity value
delta_t = 0.05
# Time interval over which the pressure drop is averaged
averaging_window = 0.025 

t_min = min(ds["t"].min() for ds in datasets.values())
t_max = max(ds["t"].max() for ds in datasets.values())
bins = np.arange(t_min, t_max + delta_t, delta_t)
bin_centers = (bins[:-1] + bins[1:]) / 2

# Compute pressure drop at each velocity step
for key, ds in datasets.items():
    mean_p, std_p = bin_and_average(ds["t"], ds["p"], bins, averaging_window)
    datasets[key]["mean"] = mean_p
    datasets[key]["std"] = std_p

#############################################################################
# Plot
#############################################################################
# Reynolds numbers for x-axis
velocity = np.minimum(0.02 * np.floor(bin_centers / 0.05 + 1), 0.3)
Re = velocity * D / nu

if args.solver and args.filter and args.drag:
    dataset_keys = list(datasets.keys())
    if args.output_suffix:
        suffix = args.output_suffix
    else:
        # Example: MB-QCM-Explicit-Implicit
        parts = [
            "-".join(s.upper() for s in args.solver),
            "-".join(f.upper() for f in args.filter),
            "-".join(d.replace("-", "").title() for d in args.drag),]
        suffix = "_".join(parts)
    dataset_keys_list = [(dataset_keys, suffix)]
else:
    dataset_keys_list = [
        (["MB PCM Explicit", "MB PCM Semi-Implicit", "MB PCM Implicit",
          "MB QCM Explicit", "MB QCM Semi-Implicit", "MB QCM Implicit"], "MB"),
        (["MF QCM Explicit", "MF QCM Semi-Implicit", "MF QCM Implicit",
          "MB QCM Explicit", "MB QCM Semi-Implicit", "MB QCM Implicit"], "MB-MF")
    ]

for keys, sub in dataset_keys_list:
    sol=plot_pressure(datasets, Re, keys, sub, args.validate)

if args.validate:
    final_solution = [Re.reshape(-1,1), sol]
    headers_list = ["Re"] + list(datasets.keys())
    
    final_solution = np.hstack(final_solution)
    file_suffix = args.output_suffix if args.output_suffix else "mf-qcm-si"
    filename = f"solution-pressure-drop-Re-{file_suffix}.dat"
    
    np.savetxt(filename, final_solution, header=" ".join(headers_list))