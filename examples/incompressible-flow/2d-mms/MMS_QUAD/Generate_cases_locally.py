"""
Summary: Script to generate different cases (as folders) of the same problem 
with different levels of refinement and polynomial degrees for the velocity and the pressure.
"""

import jinja2
import os
import numpy as np
import shutil

PATH = os.getcwd()

# User input
CASE_PREFIX = 'MMS_2D_steady_'
PRM_FILE = 'MMS_2D_steady.prm'
n_poly_deg = 3

poly_deg_1=1
poly_deg_2=3

poly_deg = np.linspace(poly_deg_1, poly_deg_2, n_poly_deg)

# Create Jinja template
templateLoader = jinja2.FileSystemLoader(searchpath=PATH)
templateEnv = jinja2.Environment(loader=templateLoader)
template = templateEnv.get_template(PRM_FILE)

# Generation of different cases

for poly_deg_u in poly_deg:
    for poly_deg_p in  range(1, int(poly_deg_u)+1):

        case_folder_name = f'{CASE_PREFIX}' + 'degu_' + f'{poly_deg_u:.0f}' + '_degp_' + f'{poly_deg_p:.0f}'


        if os.path.exists(case_folder_name) and os.path.isdir(case_folder_name):
            shutil.rmtree(case_folder_name)
    
        parameters = template.render(Poly_deg_u = int(poly_deg_u), Poly_deg_p = int(poly_deg_p))
        # Create the folder of the case and put the prm template in it
        case_path = f'{PATH}/{case_folder_name}'
        os.mkdir(case_path)
        shutil.copy(f'{PATH}/{PRM_FILE}', f'{case_path}/{PRM_FILE}')

        # Write a unique prm file with the prm template being updated
        with open(f'{case_path}/{PRM_FILE}', 'w') as f:
            f.write(parameters)

        print(f'{case_folder_name} generated successfully!')
