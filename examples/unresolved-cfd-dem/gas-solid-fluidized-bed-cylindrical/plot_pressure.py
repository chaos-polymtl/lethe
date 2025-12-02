# SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

import numpy as np
import pyvista as pv
import pandas as pd

import matplotlib.pyplot as plt
from cycler import cycler


colors=['#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02']

plt.rcParams['axes.prop_cycle'] = cycler(color=colors)

t,p,p2 = np.loadtxt("output/pressure_drop.dat",unpack=True,skiprows=1)
tt,pp,pp2 = np.loadtxt("output_modelB/pressure_drop.dat",unpack=True,skiprows=1)

D= 0.02
dp = 5e-4
nu = 1e-5
g = 9.81
rho_p = 1000
rho_f = 1
dp = 5e-4
mu = 1e-5
Ar =  g * rho_f * (rho_p-rho_f)*dp**3 / mu**2
print("Ar: ", Ar)
Re_mf_WY = (33.7**2+0.0408*Ar)**0.5 - 33.7
U_mf_WY = Re_mf_WY * mu / dp 
Re_mf_WY = U_mf_WY * D/nu
print("Wen-Yu")
print("  Re_p: ",Re_mf_WY)
print("  U_mf: ",U_mf_WY)
Re_mf_N = (19.29**2+0.0276*Ar)**0.5 - 19.29
U_mf_N = Re_mf_N * mu / dp
Re_mf_N= U_mf_N * D/nu
print("Noda")
print("  Re_p: ",Re_mf_N)
print("  U_mf: ",U_mf_N)


Re = 0.2*t * D / nu
Repp=2*tt *D/nu
exp_p=408.34
plt.plot(Re,p)
plt.plot(Repp,pp)
plt.plot([0,max(Re)],[exp_p,exp_p],'k--',label="Analytical pressure drop")
plt.plot([Re_mf_WY, Re_mf_WY],[0,450],':',label="Wen-Yu")
plt.plot([Re_mf_N, Re_mf_N],[0,450],'-.',label="Noda")
plt.xlabel("Re")
plt.ylabel("$\\Delta p$")
plt.legend()
plt.show()
