#############################################################################
"""
Code for the analytical solution of the capillary wave from Prosperetti [1].
This code solves numerically equation 22.
All units are in SI.
"""
#############################################################################

#############################################################################
'''Importing Libraries'''
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mpmath as mpm
#############################################################################

#############################################################################
# Run script: python3 capillary-wave-prosperetti-solution.py path_to_exported_analytical_data_csv_file

# Check the number of input arguments
if len(sys.argv)!=2:
    print("********************************************************************************\n"
          "Incorrect number of arguments\n"
          "Run script with: \n"
          "\t python3 capillary-wave-prosperetti-solution.py path_to_exported_analytical_data_file\n"
          "********************************************************************************\n")
    exit(1)

#############################################################################
#----------------------------------
# Physical properties
#----------------------------------

# Upper fluid properties (fluid 0)
rho_u = 1         # density
nu_u = 5e-6       # kinematic viscosity
mu_u = nu_u*rho_u # dynamic viscosity

# Lower fluid properties (fluid 1)
rho_l = 1         # density
nu_l = 5e-6       # kinematic viscosity
mu_l = nu_l*rho_l # dynamic viscosity

# Surface tension coefficient
sigma = 0.01

#----------------------------------
# Single wave dimensions
#----------------------------------

# Wavelength
l = 1e-4

# Wavenumber
k = 2*np.pi/l

# Inviscid initial angular frequency squared
omega0Sq = sigma*k**3/(rho_u+rho_l)

# Initial amplitude of the wave
a0 = 0.01*l

# Initial velocity at the center point of the wave
u0 = 0

# Characteristic timescale
tau = 1/np.sqrt(omega0Sq)

#----------------------------------
# Solution parameters
#----------------------------------

# End time
t_end = 50*tau

# Time-step
dt = tau/50

# Time array
time_array = np.arange(dt, t_end, dt)

# Output filename
filename = sys.argv[1]

#############################################################################
#----------------------------------
# Compute equation 22
#----------------------------------

lambda_u = lambda s: mpm.sqrt(k*k + s/nu_u)
lambda_l = lambda s: mpm.sqrt(k*k + s/nu_l)

Lambda = lambda s: (4*k*( - rho_l*rho_u*s
                          + k*(mu_u-mu_l)*(rho_u*(k-lambda_l(s)) - rho_l*(k-lambda_u(s)))
                          + k**2*(mu_l-mu_u)**2*(k-lambda_l(s))*(k-lambda_u(s))/s)
                      /((rho_l+rho_u)*(rho_l*(k-lambda_u(s)) + rho_u*(k-lambda_l(s)))))

A = lambda s: 1/s * (a0 + (u0*s - omega0Sq*a0) / (s**2 + Lambda(s) * s + omega0Sq))

# List of analytical amplitudes
a = []
for t in time_array:
    a.append(mpm.invlaptalbot(A,t, degree=100)) # Depending on the solution, you might need to increase the degree

relative_amplitude = [h/a0 for h in a]

#############################################################################
#----------------------------------
# Write csv file with values
#----------------------------------
lethe_df = pd.DataFrame({'t': time_array, 'a/a0': relative_amplitude})
lethe_df.to_csv(filename, index=False)

#----------------------------------
# Plot solution for visualization
#----------------------------------
fig0 = plt.figure(figsize=(9.5, 6))
ax0 = fig0.add_subplot(111)
plt.plot(time_array, relative_amplitude, "--r", linewidth=2, label=f'Prosperetti (1981) Eq22')
plt.xlabel('Time')
plt.ylabel('Relative amplitude')
plt.ticklabel_format(axis='both', style='sci', scilimits=(4,-9))
plt.legend(loc="best")
plt.tight_layout()
# plt.show() # Uncomment line to see the analytical solution

#############################################################################
#----------------------------------
# References
#----------------------------------

# [1] A. Prosperetti, “Motion of two superposed viscous fluids,” Phys. Fluids, vol. 24, no. 7, pp. 1217–1223, Jul. 1981, doi: 10.1063/1.863522.