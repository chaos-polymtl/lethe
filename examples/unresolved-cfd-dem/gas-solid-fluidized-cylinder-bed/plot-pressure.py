# SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#############################################################################
# Import Libraries
#############################################################################
import numpy as np
import pyvista as pv
import pandas as pd

import matplotlib.pyplot as plt
from cycler import cycler

#############################################################################
# Steup plot parameters
#############################################################################
colors=['#1b9e77','#d95f02','#7570b3']
markers = ['o', 's']

plt.rcParams['axes.prop_cycle'] = cycler(color=colors)

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
def plot_pressure(Re, dataset_keys, labels, subscript):
    
    for i, (key, label) in enumerate(zip(dataset_keys, labels)):
        marker = markers[(i // 3) % len(markers)]
        plt.errorbar(Re, datasets[key]["mean"], yerr=datasets[key]["std"], fmt=marker+'-',  markerfacecolor='none', markersize=8, capsize=3, label=label)
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
    plt.savefig(f"Figure_{subscript}.png", dpi=600, bbox_inches='tight')
    plt.show()


#############################################################################
# Load data
#############################################################################
data_files = {
    "MB_PCM_E": "output_mb_pcm_explicit/pressure_drop.dat",
    "MB_PCM_SI": "output_mb_pcm_semi_implicit/pressure_drop.dat",
    "MB_PCM_I": "output_mb_pcm_implicit/pressure_drop.dat",
    "MB_QCM_E": "output_mb_qcm_explicit/pressure_drop.dat",
    "MB_QCM_SI": "output_mb_qcm_semi_implicit/pressure_drop.dat",
    "MB_QCM_I": "output_mb_qcm_implicit/pressure_drop.dat",
    "MF_QCM_E": "output_mf_explicit/pressure_drop.dat",
    "MF_QCM_SI": "output_mf_semi_implicit/pressure_drop.dat",
    "MF_QCM_I": "output_mf_implicit/pressure_drop.dat"
}

datasets = {}
for key, path in data_files.items():
    t, p, _ = np.loadtxt(path, unpack=True, skiprows=1)
    datasets[key] = {"t": t, "p": p}

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

plot_pressure(
    Re,
    ["MB_PCM_E", "MB_PCM_SI", "MB_PCM_I", "MB_QCM_E", "MB_QCM_SI", "MB_QCM_I"],
    ["MB PCM Explicit", "MB PCM Semi-Implicit", "MB PCM Implicit", "MB QCM Explicit", "MB QCM Semi-Implicit", "MB QCM Implicit"],
    "MB"
)

plot_pressure(
    Re,
    ["MF_QCM_E", "MF_QCM_SI", "MF_QCM_I", "MB_QCM_E", "MB_QCM_SI", "MB_QCM_I"],
    ["MF QCM Explicit", "MF QCM Semi-Implicit", "MF QCM Implicit", "MB QCM Explicit", "MB QCM Semi-Implicit", "MB QCM Implicit"],
    "MB-MF"
)
