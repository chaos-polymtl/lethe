# SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

import numpy as np
import pyvista as pv
import pandas as pd

import matplotlib.pyplot as plt
from cycler import cycler


colors=['#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02']

plt.rcParams['axes.prop_cycle'] = cycler(color=colors)

t,p,p2 = np.loadtxt("output_modelB_test/pressure_drop.dat",unpack=True,skiprows=1)
tA,pA,p2A = np.loadtxt("output_modelA/pressure_drop.dat",unpack=True,skiprows=1)
# t,p,p2 = np.loadtxt("output_ramp/pressure_drop.dat",unpack=True,skiprows=1)

D= 0.02
# Bed cross sectional area
Ab = np.pi*(D/2)**2
dp = 5e-4
Vp = 4/3*np.pi*(dp/2)**3
Np = 200000
nu = 1e-5
g = 9.81
rho_p = 1000
rho_f = 1
mu = nu*rho_f
Ar =  g * rho_f * (rho_p-rho_f)*dp**3 / mu**2
print("Ar: ", Ar)
Re_mf_WY = (33.7**2+0.0408*Ar)**0.5 - 33.7
U_mf_WY = Re_mf_WY * nu / dp 
Re_mf_WY_inlet = U_mf_WY * D/nu
print("Wen-Yu")
print("  Re_p: ",Re_mf_WY)
print("  U_mf: ",U_mf_WY)
Re_mf_N = (19.29**2+0.0276*Ar)**0.5 - 19.29
U_mf_N = Re_mf_N * nu / dp
Re_mf_N_inlet= U_mf_N * D/nu
print("Noda")
print("  Re_p: ",Re_mf_N)
print("  U_mf: ",U_mf_N)

# Analytical pressure drop
delta_p_analytical = Np*Vp*(rho_p-rho_f)*g / Ab
print("Analytical pressure drop: ", delta_p_analytical)

# Time interval for each velocity
delta_t = 0.05
# delta_t = 0.001

exp_p=408.34

t_min = t.min()
t_max = t.max()

# Bin edges from t_min to t_max with step delta_t
bins = np.arange(t_min, t_max + delta_t, delta_t)
indices = np.digitize(t, bins) - 1

mean_p = []
std_p = []
bin_centers = []

for i in range(len(bins)-1):
    mask = indices == i
    if np.any(mask):  # check if there are points in this bin
        mean_p.append(p[mask].mean())
        std_p.append(p[mask].std())
        bin_centers.append((bins[i] + bins[i+1])/2)

bin_centers = np.array(bin_centers)
print("Bin centers:", bin_centers)
mean_p = np.array(mean_p)
std_p = np.array(std_p)
print("Mean pressures:", mean_p)

mean_p_A = []
std_p_A = []

for i in range(len(bins)-1):
    mask = indices == i
    if np.any(mask):  # check if there are points in this bin
        mean_p_A.append(pA[mask].mean())
        std_p_A.append(pA[mask].std())

bin_centers = np.array(bin_centers)
print("Bin centers:", bin_centers)
mean_p_A = np.array(mean_p_A)
std_p_A = np.array(std_p_A)
print("Mean pressures:", mean_p_A)

velocity = np.minimum(0.02 * np.floor(bin_centers / 0.05 + 1), 0.2)
# velocity = np.minimum(0.2 * bin_centers, 0.2)
Re = velocity * D / nu
print("Reynolds numbers:", Re)

plt.plot(Re, mean_p, label="Model B")
# plt.plot(Re, mean_p_A, label="Model A")
plt.plot([0,max(Re)],[delta_p_analytical,delta_p_analytical],'k--',label="Analytical pressure drop")
plt.plot([Re_mf_WY_inlet, Re_mf_WY_inlet],[0,450],':',label="Wen-Yu")
plt.plot([Re_mf_N_inlet, Re_mf_N_inlet],[0,450],'-.',label="Noda")
plt.xlabel("Re")
plt.ylabel("$\\Delta p$")
plt.legend()
plt.show()