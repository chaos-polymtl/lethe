#Post-processing tool to track mixing in small rotating drum

# importing sys
import sys

#Append path to the module to sys path list to be able to import it
path_to_module = '../../../contrib/postprocessing/' #Folder where the module is
sys.path.append(path_to_module)

print(sys.path)

#Import tools from module
from lethe_pyvista_tools import *


#Create object named "particles"
particles = lethe_pyvista_tools(case_path = ".", prm_file_name = "small-rotating-drum-dem.prm")

#Read files and store data to object named "particles"
#The packing of particles is finished after 40th time step
#So, it is important to set first to 40
particles.read_lethe_to_pyvista(pvd_name = "out.pvd", first = 40)

condition = "(y**2 + z**2)**(1/2) > 0.025"

particles.array_modifier(new_array_name = "particle_color", condition = condition, array_values = 1, restart_array = True, write_new_vtu = False)

print(particles.df_0.array_names)

condition = "(y**2 + z**2)**(1/2) > 0.04"

particles.array_modifier(new_array_name = "particle_color", condition = condition, array_values = 2, restart_array = False, write_new_vtu = True)