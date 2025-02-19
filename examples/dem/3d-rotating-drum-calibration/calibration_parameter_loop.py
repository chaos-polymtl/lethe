
import numpy as np
import jinja2
from datetime import datetime
import os
import shutil

# Loops
friction_coefficient = np.linspace(0.2,0.5,5)
rolling_friction_coefficient =  np.linspace(0,0.2,5)
restitution_coeff =np.linspace(0.8,1,3)

BASE_PREFIX = "rotating-drum"
folder = "calibration_"

for i in friction_coefficient:
    for j in rolling_friction_coefficient:
        for k in restitution_coeff:
            CASE_PREFIX = f"{int(100 * i):02d}_{int(100 * j):02d}_{int(100 * k):02d}"

            # Define the directory path based on CASE_PREFIX
            directory_path = "./prm/" + CASE_PREFIX

            # Check if the directory already exists. If so, remove it.
            if os.path.exists(directory_path):
                shutil.rmtree(directory_path)


            # Create the directory (which is now guaranteed not to exist)
            os.makedirs(directory_path)

            # Template processing and file writing
            PRM_FILE = BASE_PREFIX + ".prm"
            PATH = os.getcwd()
            templateLoader = jinja2.FileSystemLoader(searchpath=PATH)
            templateEnv = jinja2.Environment(loader=templateLoader)
            template = templateEnv.get_template(PRM_FILE)

            # Replace the values for the spreading prm
            output_text = template.render(Friction_coefficient=i,
                                          Rolling_friction_coefficient=j,
                                          Restitution_coefficient=k)

            # New prm name
            prm_file_name = BASE_PREFIX + ".prm"
            output_file_path = os.path.join(directory_path, prm_file_name)

            # Write
            with open(output_file_path, 'w') as f:
                f.write(output_text)


            print(f"{prm_file_name} has been written in {directory_path}")

shutil.copy("./particles.input","./prm")

print(f"Job is done!\n")
