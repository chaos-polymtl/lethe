#Post-processing tool to track mixing in small rotating drum

# import sys
import sys

# Path to the module
path_to_module = '../../../contrib/postprocessing/'
sys.path.append(path_to_module)

# Import tools from module
from lethe_pyvista_tools import *
import matplotlib.pyplot as plt


# Create object named "particles"
particles = lethe_pyvista_tools(".", "small-rotating-drum-dem.prm", "out.pvd")


# State condition for particle_color array creation
condition = "(y**2 + z**2)**(1/2) > 0.025"

# Create array named particle_color and give the value of 1 to all particles
# within the given condition
# The packing of particles is finished after 40th output time step
# So, it is important to set reference_time_step to 40
particles.modify_array(array_name = "particle_color", condition = condition, array_values = 1, restart_array = True, reference_time_step = 40)


# State a second condition to modify particle_color
condition = "(y**2 + z**2)**(1/2) > 0.04"

# Modify the array to give the value of 2 to all particles within the new
# condition
particles.modify_array(array_name = "particle_color", condition = condition, array_values = 2, restart_array = False, reference_time_step = 40)


# VISUALIZE DATA USING PYVISTA:

# Create a sphere with diameter 1
sphere = pv.Sphere(theta_resolution=50, phi_resolution=50)

# Use sphere as basis to create sheric representation of particles
particle_glyph = particles.get_df(40).glyph(scale='diameter', geom=sphere)

# Create a plotter object
plotter = pv.Plotter()
# Add particles to object and color them by particle_color array
plotter.add_mesh(particle_glyph, scalars = "particle_color")
# Open plot window
plotter.show()


# MIXING INDEX:

# To measure the mixing index using the nearest neighbors technique by
# Godlieb et al. (2007), the domain must be split in half. Thus:

# Find center of mass radial position
# To do so, first we will need the cylindric coods of each particle
# Note that the radial components of our cylinder are y and z.
particles.get_cylindrical_coords(radial_components = "yz")

# Since all particles are of the same size (radius = coord 0)
r_center_mass = np.mean(particles.get_df(40)["points_cyl"][:, 0])

# Split domain in half (restarting array)
condition = f"(y**2 + z**2)**(1/2) > {r_center_mass}"
particles.modify_array(array_name = "particle_color", condition = condition, array_values = 1, restart_array = True, reference_time_step = 40)

# Get list of nearest neighbors of each particle
particles.get_nearest_neighbors(return_id = True, n_neighbors = 15)

# Print position of nearest neighbors:
print(f"\nPosition of nearest neighbor of particle 2 at time-step 5 = {particles.get_df(5).points[particles.get_df(5)['neighbors'][2][0]]}\n")

# Calculate mixing index using nearest neighbors technique by
# Godlieb et al. (2007)
particles.mixing_index_nearest_neighbors(reference_array = "particle_color", n_neighbors = 15, mixing_index_array_name = "mixing_index_NNM")

# Store mixing index in a variable to compare with Doucet et al. (2008)
particles.mixing_index_nnm = particles.mixing_index

# Calculate mixing index using Doucet et al. (2008)

# Since the problem is cylindrical, we will set use_cyl = True
# Also, since Doucet does not use the color of particles, we need
# to set reference_time_step = 40
# Additionally, by default, Doucet mixing index decreases with mixing,
# but, to compare with NNM, we need to use and increasing index.
# Finally, to avoid values above 1, we normalize it using the mixing index
# of the 40th time-step seting normalize = True. 
particles.mixing_index_doucet(reference_time_step = 40, use_cyl = True, increasing_index = True, normalize = True)

# Store the second index:
particles.mixing_index_doucet = particles.mixing_index

# Save plot with mixing index as a function of time
# Note that the mixing index will only make sense for
# a time-steps greater than the one colored.
# In this case, we will plot only from the 40th time-step on.
plt.plot(particles.time_list[40:], particles.mixing_index_nnm[40:], '-b', label = "NNM")
plt.plot(particles.time_list[40:], particles.mixing_index_doucet[40:], '--k', label = "Doucet")
plt.plot(particles.time_list[40:], np.repeat(1, len(particles.time_list[40:])), ':r')
plt.xlabel("Time [s]")
plt.ylabel("Mixing index [-]")
plt.xlim(particles.time_list[40], particles.time_list[-1])
plt.ylim(0, 1.1)
plt.legend()
plt.savefig("./mixing_index.png")
plt.close()

