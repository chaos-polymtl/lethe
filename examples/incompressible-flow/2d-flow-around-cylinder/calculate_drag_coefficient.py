# SPDX-FileCopyrightText: Copyright (c) 2022-2023 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

import numpy as np

# Get the cell number and the force from the dat file
cells = np.loadtxt('force.00.dat',skiprows=1,usecols=(0))
force = np.loadtxt('force.00.dat',skiprows=1,usecols=(1))

# Calculate the drag coefficient
U_infty=1
D=1
Cd = 2*force/(U_infty**2*D)
print("Cells ", " Cd")
for i in range(0,len(Cd)):
  print(cells[i]," ", Cd[i])