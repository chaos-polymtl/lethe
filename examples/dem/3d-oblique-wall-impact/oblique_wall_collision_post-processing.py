# To use the lethe_pyvista_tools, you need to have python 3
# installed in your computer

# The modules necessary to run lethe pyvista tools are:
# Pandas: pip install pandas
# PyVista: pip install pyvista

# To use lethe_pyvista_tools, add /contrib/postprocessing to your PythonPATH 
# adding 
# ``export PYTHONPATH="${PYTHONPATH}:$LETHE_PATH/contrib/postprocessing"``
# to your ~/.bashrc, or append the path to the /contrib/postprocessing
# folder in Lethe to sys

import sys
import matplotlib.pyplot as plt
import numpy as np
# Path to the module
path_to_module = '$LETHE_PATH/contrib/postprocessing/'
sys.path.append(path_to_module)

# or even put the "lethe_pyvista_tools.py" file inside
# the same directory as your python script and procceed as follows

# This line imports all lethe_pyvista_tools functionalities
from lethe_pyvista_tools import *

# Proprieties
r_p = 0.1
h_init = 0.5
j = np.linspace(0,12,13)


# list of angles to read
n=34
theta = np.linspace(1,70,n) 
theta_out = np.zeros(n)
omega_out = np.zeros(n)
restitution_out = np.zeros(n)
vy_out = np.zeros(n)


# Experimental data from A.H. Kharaz, D.A. Gorham, and A.D. Salman. An experimental study of the elastic rebound of spheres. Powder Technology, 120(3):281 – 291, 2001. URL: http://www.sciencedirect.com/science/article/pii/S0032591001002832, doi:https://doi.org/10.1016/S0032-5910(01)00283-2.
theta_restitution_exp = [1.84, 3.703, 9.496, 20.012, 30.127, 39.382, 49.694, 60.02]
restitution_exp = [0.79344, 0.74533, 0.67257, 0.59691, 0.68756, 0.77148, 0.84761,0.90142]

theta_omega_exp = [6.143, 10.991, 21.251, 31.157, 39.143, 50.348, 60.328]
omega_exp = [136.18,265.69,577.65,612.88,595.49, 453.21, 382.09]

theta_rebound_exp =[0.016, 2.15, 4.096, 9.826, 20.567, 30.29, 39.535, 49.747, 60.208]
theta_out_exp = [0.074, 1.733, 2.962, 6.518, 12.208, 21.67, 32.613, 45.125, 57.474]

# MFIX results from https://mfix.netl.doe.gov/doc/vvuq-manual/main/html/dem/dem-05.html
theta_restitution_mfix = [0,1.602,3.203,4.522,5.297,6.12, 6.382,6.926,7.139,8.49,8.984,9.775,10.831,12.183,13.751,16.425,19.777,22.932,25.575,28.02,30.416,39.107,49.977,60.747]
restitution_mfix = [1.000,1,0.965,0.94,0.91,0.869,0.846,0.823,0.798,0.7507,0.72,0.70475,0.68,0.657,0.637,0.622,0.617,0.626,0.641,0.658,0.687,0.7773,0.844,0.897]






for i,x in enumerate(theta): # loop over the .prm
    e = lethe_pyvista_tools('.', f"run_oblique_collision_{int(x):02d}.prm", 'out.pvd',read_to_df=True)
    # Get the last vtu
    last_df = e.get_df(len(e.list_vtu)-1)
    vy = last_df['velocity'][0, 1]
    vy_out[i] = vy
    vz = last_df['velocity'][0, 2]
    omega = last_df['omega'][0, 0]

    theta_out[i] = np.arctan(vy/vz) /2 / np.pi * 360
    omega_out[i] = np.abs(omega)
    restitution_out[i] = min(np.abs(vy) / (3.9 * np.sin(x*2 * np.pi / 360)+1e-16),1)


plt.plot(theta,theta_out,'-',color="black",label="Lethe")
plt.plot(theta_rebound_exp,theta_out_exp,'s',color="black",label="Kharaz, Gorham, and Salman")
plt.legend()
plt.xlabel("Angle")
plt.ylabel("Rebound angle")
plt.savefig("rebound.png",dpi=300)
plt.show()

plt.plot(theta,omega_out,'-',color="black",label="Lethe")
plt.plot(theta_omega_exp,omega_exp,'s',color="black",label="Kharaz, Gorham, and Salman")
plt.legend()
plt.ylabel("Angular velocity")
plt.xlabel("Angle")
plt.savefig("omega.png",dpi=300)
plt.show()

plt.plot(theta,restitution_out,'-',color="black",label="Lethe")
plt.plot(theta_restitution_exp,restitution_exp,'s',color="black",label="Kharaz, Gorham, and Salman")
plt.plot(theta_restitution_mfix,restitution_mfix,'^',color="black",label="MFIX")

plt.legend()
plt.ylabel("Tangential coefficient of restitution")
plt.xlabel("Angle")
plt.savefig("coeff_restitution.png",dpi=300)
plt.show()