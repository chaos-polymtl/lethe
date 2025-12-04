# SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

import numpy as np
import pyvista as pv
import pandas as pd

import matplotlib.pyplot as plt
from cycler import cycler


colors=['#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02']

plt.rcParams['axes.prop_cycle'] = cycler(color=colors)

# Helper functions
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
def plot_pressure(Re, dataset_keys, labels):
    for key, label in zip(dataset_keys, labels):
        plt.errorbar(Re, datasets[key]["mean"], yerr=datasets[key]["std"], fmt='o-', capsize=3, label=label)
    plt.plot([0, max(Re)], [delta_p_analytical]*2, 'k--', label="Analytical Δp")
    plt.plot([Re_mf_WY_inlet]*2, [0, max([datasets[k]["mean"].max() for k in dataset_keys])*1.1], ':', label="Wen-Yu")
    plt.plot([Re_mf_N_inlet]*2, [0, max([datasets[k]["mean"].max() for k in dataset_keys])*1.1], '-.', label="Noda")
    plt.xlabel("Re")
    plt.ylabel("$\\Delta p$")
    plt.legend()
    plt.show()


# === Load data ===
data_files = {
    "A_MB": "output_modelA_mb/pressure_drop.dat",
    "B_MB": "output_modelB_mb/pressure_drop.dat",
    "A_MB_proj": "output_modelA_project_forces/pressure_drop.dat",
    "B_MB_proj": "output_modelB_project_forces/pressure_drop.dat",
    "A_MF": "output_modelA_mf/pressure_drop.dat",
    "A_MF_ramp_bdf2": "output_modelA_mf_ramp/pressure_drop.dat",
    "A_MF_ramp": "output_modelA_mf_ramp_bdf2/pressure_drop.dat"
}

datasets = {}
for key, path in data_files.items():
    t, p, _ = np.loadtxt(path, unpack=True, skiprows=1)
    datasets[key] = {"t": t, "p": p}
    
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

# Minimum fluidization velocity
Re_mf_WY = (33.7**2+0.0408*Ar)**0.5 - 33.7
U_mf_WY = Re_mf_WY * nu / dp 
Re_mf_WY_inlet = U_mf_WY * D/nu

Re_mf_N = (19.29**2+0.0276*Ar)**0.5 - 19.29
U_mf_N = Re_mf_N * nu / dp
Re_mf_N_inlet= U_mf_N * D/nu

print("Wen-Yu: Re_p =", Re_mf_WY, "U_mf =", U_mf_WY)
print("Noda: Re_p =", Re_mf_N, "U_mf =", U_mf_N)

# Analytical pressure drop
delta_p_analytical = Np*Vp*(rho_p-rho_f)*g / Ab
print("Analytical pressure drop: ", delta_p_analytical)

# Calculating the pressure drop at each velocity step for all datasets

# Time interval for each velocity
delta_t = 0.05
averaging_window = 0.025 # Time interval over which the pressure drop is averaged

t_min = min(ds["t"].min() for ds in datasets.values())
t_max = max(ds["t"].max() for ds in datasets.values())
bins = np.arange(t_min, t_max + delta_t, delta_t)
bin_centers = (bins[:-1] + bins[1:]) / 2

# Compute pressure drop at each velocity step
for key, ds in datasets.items():
    mean_p, std_p = bin_and_average(ds["t"], ds["p"], bins, averaging_window)
    datasets[key]["mean"] = mean_p
    datasets[key]["std"] = std_p

# Reynolds numbers for x-axis
velocity = np.minimum(0.02 * np.floor(bin_centers / 0.05 + 1), 0.4)
Re = velocity * D / nu

plot_pressure(
    Re,
    ["A_MB", "A_MF", "B_MB"],
    ["Model A MB", "Model A MF", "Model B MB"]
)

plot_pressure(
    Re,
    ["A_MB_proj", "A_MF", "B_MB_proj"],
    ["Model A MB-proj", "Model A MF", "Model B MB-proj"]
)
