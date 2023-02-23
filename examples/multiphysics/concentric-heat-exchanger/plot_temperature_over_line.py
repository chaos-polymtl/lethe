"""
Postprocessing code for the concentric heat exchanger to extract the temperature profile

Author : Bruno Blais
"""

# Modules
#-------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv

#--------------------------------------------
# Main
#--------------------------------------------

# Load VTU file
vtu_file="./output/out.0001.0000.vtu"
sim = pv.read(vtu_file)
sim.set_active_scalars("temperature")


# Create begin and end point of line
a = [0., 0, 0]
b = [0., 0, 50]
# Extract all field over the line using pyvista
sampled_data_center=sim.sample_over_line(a, b, resolution=1000)

# Create begin and end point of line
a = [1, 0, 0]
b = [1, 0, 50]
sampled_data_periphery=sim.sample_over_line(a, b, resolution=1000)

# Create begin and end point of line
a = [0.5, 0, 0]
b = [0.5, 0, 50]
sampled_data_mid=sim.sample_over_line(a, b, resolution=1000)

# Extract temperature along sampled lines
z_c = sampled_data_center["Distance"]
T_c = sampled_data_center["temperature"][:]

z_p = sampled_data_periphery["Distance"]
T_p = sampled_data_periphery["temperature"][:]

z_m = sampled_data_mid["Distance"]
T_m = sampled_data_mid["temperature"][:]

plt.plot(z_c,T_c,label="Center of channel")
plt.plot(z_m,T_m,label="Half radius of channel")
plt.plot(z_p,T_p,label="Inner wall")
plt.plot([50],[25.25],'ko',label="NTU approach with correlations")

plt.xlabel("z [mm]")
plt.ylabel("Temperature [$^\circ$C]")
plt.legend()

plt.savefig("temperature_along_line.png",dpi=200)
plt.show()
