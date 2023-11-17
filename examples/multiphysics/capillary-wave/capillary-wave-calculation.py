#############################################################################
"""
Script for the capillary wave example: compute simulation end time and
capillary time-step constraint.
"""
#############################################################################

#############################################################################
'''Importing Libraries'''
import numpy as np
import sys
#############################################################################

#############################################################################
# Run script: python3 capillary-wave-calculation.py
# or
# Run script: python3 capillary-wave-calculation.py path_to_exported_file_with_calculations

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

rho_sum = rho_u + rho_l
mu_sum = mu_u+ mu_l

# Surface tension coefficient
sigma = 0.01

#----------------------------------
# Single wave dimensions
#----------------------------------
# Wavelength
l = 1e-4

# Wavenumber
k = 2*np.pi/l

# Initial amplitude of the wave
a0 = 0.01*l

# Initial angular frequency of the wave (inviscid fluid assumption)
omega0 = np.sqrt(sigma*k**3/rho_sum)

#----------------------------------
# Simulation parameters
#----------------------------------
# Desired mesh resolution
desired_dx = 0.01*l

# Initial subdivision of the domain along the x-axis
initial_subdivision = 4

# Compute desired number of refinement
n_refinement_desired = int(np.ceil(np.log(l/desired_dx/initial_subdivision)/np.log(2)))

# Imposed refinement
n_refinement = 5 # Make sure that this value corresponds to the finest refinement level of your simulation

#############################################################################
#----------------------------------
# Real resolution of the mesh
#----------------------------------
real_dx = l/(initial_subdivision*2**n_refinement)

#----------------------------------
# Ohnesorge number
#----------------------------------
Oh = mu_sum/np.sqrt(rho_sum*sigma*2*real_dx)

#-------------------------------------
# Compute end time and capillary
# time-step constraint
#-------------------------------------
t_end = 50/omega0
time_step_constraint = np.sqrt(rho_sum*real_dx**3/(2*np.pi*sigma))

#############################################################################
#-------------------------------------
# Print computation results
#-------------------------------------
if len(sys.argv)<2:
    print(f" End time: {t_end} s")
    print(f" The maximal mesh resolution should be: {desired_dx} m")
    print(f" The number of refinement should be: {n_refinement_desired}")
    print(f" The number of refinement imposed is: {n_refinement}")
    print(f" The mesh resolution is: {real_dx} m")
    print(f" The capillary constraint is: {time_step_constraint} s")
    print(f" Oh = {Oh}")
else:
    output_file = sys.argv[1]
    print(f" Computed values are saved in: {output_file}")
    with open(output_file, 'w') as file:
        t_end = 50/omega0
        print(f"End time: {t_end} s",file=file)
        print(f"The maximal mesh resolution should be: {desired_dx} m",file=file)
        print(f"The number of refinement should be: {n_refinement_desired}",file=file)
        print(f"The number of refinement imposed is: {n_refinement}",file=file)
        print(f"The mesh resolution is: {real_dx} m",file=file)
        print(f"The capillary constraint is: {time_step_constraint} s",file=file)
        print(f"Oh = {Oh}",file=file)