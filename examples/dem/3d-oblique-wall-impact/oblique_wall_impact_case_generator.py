import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy import optimize
import jinja2
import os

# Case generator
PATH = os.getcwd()
PRM_FILE = 'oblique_wall_impact_template.prm'

# System characteristics
theta = np.linspace(1,70,34)      # Restitution coefficient
velocity = 3.9 

templateLoader = jinja2.FileSystemLoader(searchpath=PATH)
templateEnv = jinja2.Environment(loader=templateLoader)
template = templateEnv.get_template(PRM_FILE)

for val in theta:

    output_text = template.render(FN=f"{int(val):02d}", vz=-velocity*np.cos(val/360*2*np.pi), vy=velocity*np.sin(val/360*2*np.pi))
    prm_file_name = f"run_oblique_impact_{int(val):02d}.prm"

    # Write the output text to the prm file
    output_file_path = os.path.join("./", prm_file_name)
    with open(output_file_path, 'w') as f:
        f.write(output_text)

