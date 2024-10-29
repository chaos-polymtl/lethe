"""
Summary: Script to generate different cases (as folders) of the same problem 
with changes in multiple parameters
"""

import jinja2
import os
import numpy as np
import shutil
import sys

PATH = os.getcwd()

# User input
CASE_PREFIX = ''
PRM_FILE = 'bubble-detachment-shear-flow.prm'
SHELL_FILE = 'launch_lethe.sh'

# Generation of data points

air_flow = [5e-7]
generic_cases = [1]
shear_rate = [100,200,300,450]
   
velocity_liquid = [0.005*shear for shear in shear_rate]
rho_w,mu_w,sigma_w=1000,1.0016e-3,0.073
L,R_0,rho_a,g,dx = 0.005,0.0005,1.23,9.81,2.443e-5
D_0 = 2*R_0
nu_w = mu_w/rho_w
velocity_air = [2*q/(np.pi*R_0*R_0) for q in air_flow]
shear_rate_air = [0.5*v_a/R_0 for v_a in velocity_air]
density = [rho_w*case for case in generic_cases]
surface_tension = [sigma_w*case for case in generic_cases]

Reynolds = []  #  rho_l*abs(U_l-U_g)*D_0/mu_l

number_of_cases = 0


# Create Jinja template
templateLoader = jinja2.FileSystemLoader(searchpath=PATH)
templateEnv = jinja2.Environment(loader=templateLoader)
template = templateEnv.get_template(PRM_FILE)

# Generation of different cases (density sweep then surface tension sweep)
if os.path.exists(f'{PATH}/shear_rate_sweep/'):
    print(f'{PATH}/shear_rate_sweep/')
    shutil.rmtree(f'{PATH}/shear_rate_sweep/')
    
    
os.mkdir(f'{PATH}/shear_rate_sweep/')

summary_sweep_path = 'summary_sweep.dat'

#Density sweep
with open(summary_sweep_path,'w') as file:
    file.write('###########################################'+ '\n')
    file.write('Shear rate sweep'+ '\n')
    file.write('###########################################'+ '\n')
    file.write('rho_l nu_l sigma v_a v_l Re capillary_dt mobility'+ '\n')
    for rho_l in density:
        for v_a in velocity_air:
            for v_l in velocity_liquid:
                number_of_cases += 1
                # Compute the Reynolds number of the case

                Re = rho_l*np.abs(v_a-v_l)*D_0/mu_w
                
                # Compute the constrained parameters according to the physical properties of the problem : capillary timestep and mobility coefficient
                capillary_dt = 2*np.sqrt(((rho_l+rho_a)*dx*dx*dx)/(4*np.pi*sigma_w))
                shear_rate_air = 0.5*v_a/R_0
                mobility = ((v_l/0.005)+shear_rate_air)*R_0*(dx*2)**2/sigma_w # Mobility suggested by Yue et al. D = S*R_0*epsilon²/sigma

                file.write(f'{rho_l:.3E} {nu_w:.3E} {sigma_w:.3E} {v_a:.3E} {v_l:.3E} {Re:.3E} {capillary_dt:.3E} {mobility:.3E}' + '\n')

                Reynolds.append(Re)

                case_folder_name = f'v_a_{v_a:.3E}_v_l_{v_l:.3E}_rho_l_{rho_l:.3E}_nu_l_{nu_w:.3E}_sigma_{sigma_w:.3E}'

                if os.path.exists(case_folder_name) and os.path.isdir(case_folder_name):
                    shutil.rmtree(case_folder_name)

                # Insert the parameters in the prm template with Jinja2 and render it
                parameters = template.render(v_a=v_a, v_l=v_l, rho_l=rho_l, nu_l=nu_w, sigma=sigma_w,capillary_dt=capillary_dt,mobility=mobility)

                # Create the folder of the case and put the prm template in it
                case_path = f'{PATH}/shear_rate_sweep/{case_folder_name}'

                os.mkdir(case_path)
                shutil.copy(f'{PATH}/{PRM_FILE}', f'{case_path}/{PRM_FILE}')

                # Write a unique prm file with the prm template being updated
                with open(f'{case_path}/{PRM_FILE}', 'w') as f:
                    f.write(parameters)
                    
                # Copy the shell script 
                shutil.copy(f'{PATH}/{SHELL_FILE}', f'{case_path}/{SHELL_FILE}') 
                
                print(f'{case_folder_name} generated successfully!')
                    


print(f'Max Re = {np.max(Reynolds):.3E}')
print(f'Min Re = {np.min(Reynolds):.3E}')

print(f'Number of cases = {number_of_cases}')
