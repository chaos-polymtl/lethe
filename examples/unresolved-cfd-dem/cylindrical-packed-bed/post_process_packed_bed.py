# SPDX-FileCopyrightText: Copyright (c) 2022-2026 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

"""
Postprocessing code for cylindrical packed bed
This code extracts the pressure along the x axis at the center of the cylinder

"""

# Modules
#-------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv
import pandas as pd

# Plot font and colors
#---------------------
font = {'weight' : 'normal',
        'size'   : 13}

plt.rc('font', **font)
colors=['#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02']

#--------------------------------------------
# Main
#--------------------------------------------
  
# Particle properties
dp=0.001
n_p=10000
rhop=2500

# Fluid properties
nu_f=0.00001
rho_f=1
mu_f=nu_f*rho_f

# Fluid inlet velocity (m/s)
U=0.2

# Column diameter
dc=0.02

# Packed bed height
# Load lethe particles last vtu file of the dem simulation
particles_vtu_file = "./output_dem/out.60000.00000.vtu"
sim_particles = pv.read(particles_vtu_file)
particles_loc = pd.DataFrame(np.copy(sim_particles.points), columns=['x', 'y','z'])
H_bed = particles_loc['x'].max() + 0.01 # The floating wall is placed at x=-0.01 

print("Packed bed height from particles: ", H_bed)

# Packed bed height
V_particles=n_p*(4/3)*np.pi*(dp/2)**3
A_bed=np.pi*(dc/2)**2
eps_packed = 1 - V_particles/(A_bed*H_bed)

print("Packed bed porosity: ", eps_packed)

# Pressure drop from the Ergun equation
deltaP_ergun= (150*(1-eps_packed)**2*mu_f*U/(dp**2*eps_packed**3) + 1.75*(1-eps_packed)*rho_f*U**2/(dp*eps_packed**3))*H_bed
print("Pressure drop from Ergun equation: ",deltaP_ergun)

# Load VTU file
vtu_file="./output/out.00001.00000.vtu"
sim = pv.read(vtu_file)
sim.set_active_scalars("pressure")

# Create begin and end point of line
a = [-0.1, 0, 0]
b = [0.1, 0, 0]

# Extract all field over the line using pyvista
sampled_data=sim.sample_over_line(a, b, resolution=1000)

# Get pressure from sampled data
x = sampled_data["Distance"]
p = sampled_data["pressure"][:]



# Get void fraction from sampled data
sim.set_active_scalars("void_fraction")
sampled_data=sim.sample_over_line(a, b, resolution=1000)
x_2 = sampled_data["Distance"]
void_fraction = sampled_data["void_fraction"][:]

# Create the figure and the first axis
fig, ax1 = plt.subplots()

# Plot the data for the first y-axis
ax1.plot(x,p,label="Pressure drop",color=colors[0])
ax1.plot([0,max(x)],[deltaP_ergun,deltaP_ergun],'k--',label="Ergun correlation")
ax1.set_xlabel('Position')  # Common x-axis label
ax1.set_ylabel('Pressure', color=colors[0])
ax1.set_ylim([-1,60])
ax1.tick_params(axis='y', labelcolor=colors[0])  # Set y-axis tick color
ax1.set_xticks([0,0.05,0.1,0.15,0.2]) #xticks(np.arange(min(x), max(x)+1, 1.0))


# Add a nice arrow label for the Ergun correlation
# Arrow + label
x_mid = 0.15 * max(x)

ax1.annotate(
    "Ergun correlation",
    xy=(x_mid, deltaP_ergun),          # point on the line
    xytext=(x_mid, 0.85 * deltaP_ergun),  # text position (slightly above)
    textcoords="data",
    ha="center",
    arrowprops=dict(arrowstyle="->", linewidth=1)
)

# Create a second y-axis sharing the same x-axis
ax2 = ax1.twinx()
ax2.plot(x_2,void_fraction,label="Void Fraction",color=colors[1])
ax2.set_ylabel('Void fraction', color=colors[1])
ax2.tick_params(axis='y', labelcolor=colors[1])  # Set y-axis tick color
plt.savefig("pressure_drop_void_fraction.png",dpi=200)
plt.show()


# Get u component of the velocity from sampled data
sim.set_active_vectors("velocity")
y = sampled_data["Distance"]
u = sampled_data["velocity"][:,0]

# Create the figure and the first axis
fig, ax1 = plt.subplots()

# Plot the data for the first y-axis
ax1.plot(x,u,label="Velocity",color=colors[2])
ax1.set_xlabel('Position')  # Common x-axis label
ax1.set_ylabel('Velocity', color=colors[2])
ax1.tick_params(axis='y', labelcolor=colors[2])  # Set y-axis tick color
ax1.set_xticks([0,0.05,0.1,0.15,0.2]) #xticks(np.arange(min(x), max(x)+1, 1.0))

# Create a second y-axis sharing the same x-axis
ax2 = ax1.twinx()
ax2.plot(x_2,void_fraction,label="Void Fraction",color=colors[1])
ax2.set_ylabel('Void fraction', color=colors[1])
ax2.tick_params(axis='y', labelcolor=colors[1])  # Set y-axis tick color
plt.savefig("velocity_void_fraction.png",dpi=200)
plt.show()

