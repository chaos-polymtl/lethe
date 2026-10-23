# SPDX-FileCopyrightText: Copyright (c) 2025-2025 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

"""
Postprocessing code for cylindrical fin
"""

# Modules
#-------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv
import pandas as pd
import argparse
import glob
import os

#Plot font and colors
font = {'weight' : 'normal',
        'size'   : 13}

plt.rc('font', **font)


#--------------------------------------------
# Main
#--------------------------------------------

parser = argparse.ArgumentParser(description='Arguments for the validation of the fin')
parser.add_argument("--vtu", type=str, help="VTU file", required=True)
parser.add_argument("--validate", action="store_true", help="Launches the script in validation mode. This will log the content of the graph and prevent the display of figures", default=False)
args, leftovers=parser.parse_known_args()

#--------------------------------------------
# Analytical solution
#--------------------------------------------
h = 10
k = 100
L = 0.2
R = 0.01
P = 2*np.pi*R
A = np.pi*R**2
T_inf = 20
T_b = 100
x_analytical = np.linspace(0,L,100)
m = np.sqrt(P*h/k/A)
T_analytical = T_inf + (T_b-T_inf) * np.cosh(m*(L-x_analytical))/np.cosh(m*L)

analytical_flux = k*A*m*(T_b-T_inf)*np.tanh(m*L)
print("Analytical heat flux: ", analytical_flux)

#--------------------------------------------
# Simulation Data
# Load VTU file
#--------------------------------------------
vtu_file=args.vtu
sim = pv.read(vtu_file)
sim.set_active_scalars("temperature")

# Create begin and end point of line
a = [-0.1, 0, 0]
b = [0.1, 0, 0]

# Extract all field over the line using pyvista
sampled_data=sim.sample_over_line(a, b, resolution=50)

# Get u component of the velocity from sampled data
x = sampled_data["Distance"]
T = sampled_data["temperature"][:]

plt.plot(x/L,T,'o',mfc='None',label="Lethe")
plt.plot(x_analytical/L,T_analytical,color="black",label="Analytical solution")

plt.xlabel("x/L")
plt.ylabel("T")
plt.legend()
if (not args.validate):
    plt.savefig("temperature_fin.png",dpi=200)
    plt.show()
else:
    plt.savefig("temperature_fin.pdf")

    # Store the temperature profile and the analytical one at the same
    # positions so that they can be compared between validation runs
    T_analytical_x = T_inf + (T_b-T_inf) * np.cosh(m*(L-x))/np.cosh(m*L)
    np.savetxt("solution-temperature.dat",
               np.column_stack((x/L, T, T_analytical_x)),
               header="x/L T T_analytical")

    # Store the heat flux on every boundary, which Lethe writes next to the
    # VTU file, with the analytical heat flux
    heat_flux_files = glob.glob(os.path.join(os.path.dirname(vtu_file), "heat_flux*.dat"))
    if (len(heat_flux_files) == 0):
        print("No heat flux file was found next to the VTU file")
    else:
        heat_flux = pd.read_csv(heat_flux_files[0], sep=r"\s+")
        flux_columns = [c for c in heat_flux.columns if c.startswith("total_flux_bc_")]
        np.savetxt("solution-heat-flux.dat",
                   np.array([[int(c.replace("total_flux_bc_", "")),
                              heat_flux[c].iloc[-1], analytical_flux]
                             for c in flux_columns]),
                   header="boundary_id total_flux analytical_flux")