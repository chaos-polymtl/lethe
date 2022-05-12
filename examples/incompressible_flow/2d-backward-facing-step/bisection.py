"""
Postprocessing code for 2D-backward-facing-step example
Computes the reattachment length (x_r) for Re = 100
for several meshes with a bisection algorithm

Author : Charles Le Pailleur
"""

import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv

########################################
########################################

# EXAMPLE TO COMPUTE x_r (REATTACHMENT LENGTH) FOR Re = 100 
# NOTE : THIS FILE MUST BE IN THE "2d-backward-facing-step" DIRECTORY TO 
#        WORK PROPERLY

# VARIABLES
L_in = 15                    # Inlet length
x_r_ref = 2.922              # Benchmark value for Re = 100 (Erturk 2008)
							  # (see attached .txt file for other Reynolds
							  #  x_r values)
x_inf = 2.                   # Initial guesses
x_sup = 3.
tol = 1e-12                  # Bisection stop criterion

# DATA EXTRACTION
n = 11                       # Number of VTU files to be read
file = list(range(0,n))
data = list(range(0,n))
N_elements = np.zeros(n)
for i in range(0,n):
    file[i] = ('Reynolds100-600/backward_facing_step_output.' 
               + f'{i:04d}' + '.0000.vtu')
    data[i] = pv.read(file[i])
    data[i].set_active_vectors("velocity")
    N_elements[i] = data[i].GetNumberOfElements(0)
    
# BISECTION METHOD
x_r = np.zeros(n)            # Reattachment length
for i in range(0,n):
    # Initial guesses
    x1 = L_in + x_inf
    x2 = L_in + x_sup
    j = 0
    # Bisection loop
    while abs(x2-x1)>tol and j<=20:
        xm = (x1+x2)/2
        # Profiles calculation
        a = np.array([x1, 0, 0])
        b = np.array([x1, 0.01, 0])
        profil1 = data[i].sample_over_line(a, b, resolution=1000)
        u1 = profil1["velocity"][:,0]
        a = np.array([xm, 0, 0])
        b = np.array([xm, 2, 0])
        profilm = data[i].sample_over_line(a, b, resolution=1000)
        um = profilm["velocity"][:,0]
        y = profil1["Distance"]
        # Derivatives
        dudy1 = (-u1[2] + 4*u1[1] - 3*u1[0]) / (y[2] - y[0])
        dudym = (-um[2] + 4*um[1] - 3*um[0]) / (y[2] - y[0])
        # Signs check
        if dudy1*dudym<0:
            x2 = xm
        else:
            x1 = xm
        j+=1
    x_r[i] = xm - L_in           # Distance from the step
    print('Mesh refinement %.f : %.4f' % (i,x_r[i]))
    
# CONVERGENCE GRAPH
mesh= np.linspace(0, n-1, n)
plt.figure()
plt.rc('axes', labelsize=15) 
plt.scatter(mesh[1:], x_r[1:])
plt.grid()
plt.xticks([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
plt.xlabel("Mesh refinement step")
plt.ylabel(r'$x_r$/h')
plt.title(r'Convergence of $x_r$ according to mesh refinement level')
plt.savefig('xr_convergence.png')

# RELATIVE ERROR
e = abs(x_r_ref-x_r)/x_r_ref

# ERROR ANALYSIS
plt.figure()
plt.loglog(N_elements[1:],e[1:])
plt.grid(which='minor')
plt.xlabel("Number of elements", fontsize=12)
plt.ylabel(r'$\frac{x_r - \bar{x}_r}{x_r}$')
plt.title('Evolution of the error with the number of elements')
plt.tight_layout()
plt.savefig('error_analysis.png')

