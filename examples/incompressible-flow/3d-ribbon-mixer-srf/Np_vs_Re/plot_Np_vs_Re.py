# SPDX-FileCopyrightText: Copyright (c) 2022, 2026 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

import numpy as np
import matplotlib.pyplot as plt

# Simulation Data
C,Tx,Ty,Tz = np.loadtxt("gather.dat",unpack=True)
n_dat = np.size(Tz)
D = 0.27
N = 10 / 2 / np.pi
Re = np.logspace(-1,np.log10(100),25)
Np = 2 * np.pi * np.abs(Tz) / 1 / N**2 / D**5


# Experimental data: three (Re, Np) series in a tab-separated file that uses
# decimal commas and in which some values are missing
experimental = [([], []) for _ in range(3)]
with open("experimental-data.dat", encoding="utf-8") as f:
    for line in f.read().splitlines()[1:]:
        columns = line.replace(",", ".").split("\t")
        for i in range(3):
            if 2*i+1 < len(columns) and columns[2*i].strip() and columns[2*i+1].strip():
                experimental[i][0].append(float(columns[2*i]))
                experimental[i][1].append(float(columns[2*i+1]))
(Re_1,Np_1), (Re_2,Np_2), (Re_3,Np_3) = experimental


# Plot data
plt.loglog(Re_1,Np_1,'o',mfc='None',label="Experiment 1")
plt.loglog(Re_2,Np_2,'*',mfc='None',label="Experiment 2")
plt.loglog(Re_3,Np_3,'^',mfc='None',label="Experiment 3")
plt.loglog(Re[:n_dat],Np,'-',label="Lethe",color="black")

plt.xlabel("$Re$")
plt.ylabel("$N_p$")
plt.legend()

plt.savefig("./Np_vs_Re.png", dpi=300)
plt.show()
