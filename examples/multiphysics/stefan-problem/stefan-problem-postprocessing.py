import numpy as np 
from scipy.optimize import fsolve
from scipy.special import erf
import matplotlib.pyplot as plt
import pyvista as pv

#Physical properties
Cp=1
dT=1
hL=100

# Value of the Stefan number
St=Cp*dT/hL

# Time at the end of the simulation
t=5

# Function to calculate the beta constant of the stefan problem
def stefan(beta):
   return beta * np.exp(beta**2) * erf(beta) - St / np.sqrt(np.pi)

# From there we can derive the position of the interface
beta = fsolve(stefan, [0.1],xtol=1e-12)
print("Beta value is : ", beta)
delta = 2 * beta * np.sqrt(t)

# Load VTU file
vtu_file="./output/stefan.0250.0000.vtu"
sim = pv.read(vtu_file)
sim.set_active_scalars("temperature")


# Create begin and end point of line
a = [0, 0.005, 0]
b = [1, 0.005, 0]

# Extract the temperature
sampled_data=sim.sample_over_line(a, b, resolution=1000)
x = sampled_data["Distance"]
T = sampled_data["temperature"]

# Calculate analytical solution. First part if the Stefan profile, after that a constant
T_analytical = 1. - (erf(x /np.sqrt(t)*0.5) / erf(beta[0]))
T_analytical[np.where(x>delta)]=0

# Make a plot out of it
plt.rcParams.update({'font.size': 15})
plt.plot(x,T,label="Lethe",lw=2)
plt.plot(x,T_analytical,"--",lw=2,label="Analytical solution")
plt.xlabel("x")
plt.ylabel("T")
plt.legend()
plt.savefig("lethe-stefan-comparison.png", dpi=300)
plt.show()


