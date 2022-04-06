"""
Postprocessing code for 2D-backward-facing-step example
Computes velocity distributions at inlet and outlet and
compares it with analytical solution (Poiseuille)

Author : Charles Le Pailleur
"""

import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv

########################################
########################################

# EXAMPLE TO COMPUTE VELOCITY DISTRIBUTION AT INLET AND OUTLET
# NOTE : THIS FILE MUST BE IN THE "2d-backward-facing-step" DIRECTORY TO 
#        WORK PROPERLY

# VARIABLES
L_out = 50                      # Outlet length
L_in = 15                       # Inlet length
h_in = 0.5                      # Half height of inlet
h_out = 1                       # Half height of outlet
u_moy_in = 1                    # Inflow mean velocity
u_moy_out = 0.5                 # Outflow mean velocity
plt.rc('axes', labelsize=14)
L_out = L_in + L_out

# DATA EXTRACTION
n = 11                          # Number of VTU files to be read
file = list(range(0,n))
data = list(range(0,n))
N = np.zeros(n)
for i in range(0,n):
    file[i] = ('Reynolds100-600/backward_facing_step_output.' 
               + f'{i:04d}' + '.0000.vtu')
    data[i] = pv.read(file[i])
    data[i].set_active_vectors("velocity")

###########################
# Poiseuille Inlet
###########################

# EXACT SOLUTION
k_in = -3/2*u_moy_in/h_in**2
y_in_an = np.linspace(-h_in, h_in, 1001)
u_in_an = k_in*(y_in_an**2-h_in**2)

# WITH LETHE
a = np.array([L_in, 1, 0])
b = np.array([L_in, 2, 0])
u_inlet_moy = np.zeros(n)
u_inlet_max = np.zeros(n)
mesh = np.linspace(0, n-1, n)
plt.figure()
for i in range(0,n):
    # Data
    data_inlet = data[i].sample_over_line(a, b, resolution=1000)
    y_inlet = data_inlet["Distance"]
    u_inlet = data_inlet["velocity"][:,0]
    # Graph
    plt.plot(u_inlet ,y_inlet, label="Mesh "+str(i))
    # Mean and max velocities
    u_inlet_moy[i] = np.mean(u_inlet)
    u_inlet_max[i] = max(u_inlet)
    
# Graph setup
plt.plot(u_in_an, y_inlet, label="Poiseuille", linestyle='--',
         linewidth=2.0, color='k')
plt.xlabel(r'$u_{in}$')
plt.ylabel("y/h")
plt.legend(fontsize=8)
plt.axis([0, 2.5, 0, 1])
plt.title(r'Velocity distribution $u(y)$ at the step')
plt.savefig('velocity_step.png')

####################
# Poiseuille Outlet
####################

# EXACT SOLUTION
k_out = -3/2*u_moy_out/h_out**2
y_out_an = np.linspace(-h_out, h_out, 1001)
u_out_an = k_out*(y_out_an**2-h_out**2)

# LETHE
a = np.array([L_out, 0, 0])
b = np.array([L_out, 2, 0])
u_outlet_moy = np.zeros(n)
u_outlet_max = np.zeros(n)
plt.figure()
for i in range(1,n):
    # Data
    data_outlet = data[i].sample_over_line(a, b, resolution=1000)
    y_outlet = data_outlet["Distance"]
    u_outlet = data_outlet["velocity"][:,0]
    # Graph
    plt.plot(u_outlet, y_outlet, label="Mesh "+str(i))
    # Mean and max velocities
    u_outlet_moy[i] = np.mean(u_outlet)
    u_outlet_max[i] = max(u_outlet)
    
# Graph setup
plt.plot(u_out_an, y_outlet, label="Poiseuille", linestyle='--',
         linewidth=2.0, color='k')
plt.xlabel(r'$u_{out}$')
plt.ylabel("y/h")
plt.legend(fontsize=8)
plt.axis([0,1.5,0,2])
plt.title(r'Velocity distribution $u(y)$ at the outlet')
plt.savefig('velocity_outlet.png')