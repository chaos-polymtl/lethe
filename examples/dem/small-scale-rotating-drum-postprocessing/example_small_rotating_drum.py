#Post-processing tool to track mixing in small rotating drum

# import sys
import sys

#Path to the module
path_to_module = '../../../contrib/postprocessing/'
sys.path.append(path_to_module)

#Import tools from module
from lethe_pyvista_tools import *


#Create object named "particles"
particles = lethe_pyvista_tools(".", "small-rotating-drum-dem.prm")

#Read files and store data to object named "particles"
particles.read_lethe_to_pyvista("out.pvd")


#State condition for particle_color array creation
condition = "(y**2 + z**2)**(1/2) > 0.025"

#Create array named particle_color and give the value of 1 to all particles within the given condition
#The packing of particles is finished after 40th time step
#So, it is important to set reference_time_step to 40
particles.array_modifier(array_name = "particle_color", condition = condition, array_values = 1, restart_array = True, reference_time_step = 40)


#State a second condition to modify particle_color
condition = "(y**2 + z**2)**(1/2) > 0.04"

#Modify the array to give the value of 2 to all particles withing the new condition
particles.array_modifier(array_name = "particle_color", condition = condition, array_values = 2, restart_array = False, reference_time_step = 40)


#VISUALIZE DATA USING PYVISTA:

#Create a sphere with diameter 1
sphere = pv.Sphere(theta_resolution=50, phi_resolution=50)

#Use sphere as basis to create sheric representation of particles
particle_glyph = particles.df[0].glyph(scale='Diameter', geom=sphere)

#Create a plotter object
plt = pv.Plotter()
#Add particles to object and color them by particle_color array
plt.add_mesh(particle_glyph, scalars = "particle_color")
#Open plot window
plt.show()

#Save results to be able to open them on ParaView
particles.write_vtu(prefix = "mod_")