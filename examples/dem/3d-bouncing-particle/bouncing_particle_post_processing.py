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


# .prm file to read
list1 = np.array(["bouncing_particle_05",
                  "bouncing_particle_06",
                  "bouncing_particle_07",
                  "bouncing_particle_08",
                  "bouncing_particle_09",
                  "bouncing_particle_10"])

color_index= 0
color_list = np.array(['r', 'b','c' ,'y' ,'g' ,'m'])
order = np.empty(len(list1)*2 , int)

for x in list1: # loop over the .prm
    e = lethe_pyvista_tools('.', x, 'out.pvd')

    # plot marker
    line = "-x" + color_list[color_index]
    square = "o" + color_list[color_index]

    position_z = np.array([])  # Initialize an empty array to store the z position at each time step
    vitesse_z = np.array([])  # Initialize an empty array to store the speed at each time step
    top_bounce = np.array([0.5])  # Initialize an empty array to store the top of each bounce
    time = e.time_list  # Initialize an empty array to store the time step

    for i in in range(len(e.list_vtu)): # Extract position and speed at every time step

        df = e.get_df(i)
        position_z = np.append(position_z,df.points[0][2])
        vitesse_z = np.append(vitesse_z,df['velocity'][0, 2])

    k = np.array([0])

    R_c = []
    for i in range(13):
        R_c = np.append(R_c, (h_init-r_p)*e.prm_dict["restitution coefficient particles"]**(2*i) + r_p )

    plt.plot(j, R_c, line, ms =12, mew = 2, label = "Analytical, e = " + str(e.prm_dict["restitution coefficient particles"]))
    order[color_index] = int(2 * (color_index))


    n = 1
    for i in range(len(vitesse_z)-1):
        if (vitesse_z[i]*vitesse_z[i+1] < 0) and (vitesse_z[i+1]<0):  # If the speed is around zero and derivative is negative
            top_bounce = np.append(top_bounce, position_z[i])
            k = np.append(k, n)
            n = n+1

    plt.plot(k, top_bounce, square, ms =6 ,mew = 1 , mec = "k", label = "Numerical, e = " + str(e.prm_dict["restitution coefficient particles"]))
    order[color_index + len(list1)] = int(order[color_index] + 1)
    color_index += 1


R_e = e.prm_dict["diameter"] / 2
Y_e =( ((1 - e.prm_dict["poisson ratio particles"]**2. ) / e.prm_dict["young modulus particles"]) + ( (1 -  e.prm_dict["poisson ratio wall"]**2 ) / e.prm_dict["young modulus wall"] ) )**(-1)
m_e = e.prm_dict["density particles"] * 4 * np.pi * R_e**3 / 3
V = 1.

k_n = np.round((16/15) * np.sqrt(R_e) * Y_e * ( (15 * m_e * V * V) / (16 * np.sqrt(R_e) * Y_e))**(0.2))

handles, labels = plt.gca().get_legend_handles_labels()
plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order],ncol=2,bbox_to_anchor=(0.3, 0.5), fontsize = 9)

plt.grid()
plt.xlabel("Bounce number")
plt.ylabel("Height of the bounce")
plt.title( r"$k_n$ = " + f"{k_n:.0e}" + r" $N \cdot m^{-1}$" )
plt.xlim(0, 12.5)
plt.show()