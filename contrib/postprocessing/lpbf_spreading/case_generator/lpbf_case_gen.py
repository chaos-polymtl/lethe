import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy import optimize
import jinja2
import argparse
import os

# Simulation parameters
coater_speed = 0.100            # Speed of the blades
number_of_layers = 20           # Number of layer. We do +2 in the code for layer -1 and 0.
diameters_list = 50E-6          # List of all the diameters we want to use
#diameter_probabilities = 1.
delta_b_p = 100E-6              # Distance between the tip of the blade and the transfert plate
delta_miss = 0.                 # Miss-match between the build-plate and the seperators
delta_o = 500E-6                # Thickness of the first layer
delta_n = 100E-6                # Thickness of the following layers
first_layer_extrusion = 2000E-6 # Extrusion of the first layer
other_layer_extrusion = 600E-6  # Extrusion of the following layers
GAP = 500E-6                    # Void around the build plate
lenght_domain = 0.0276          # Lenght of the domain. Is use to compute the time it takes the blades to move
alloy            = "ALUM"       # Input for case generator
first_starting_time = 0.3                        # --> HARD CODED IN THE POST_PROCESSING. This could change depending on the number of particles being used.
plates_speed = 0.005                             # Y speed of both plate
plate_displacement_time = delta_n/plates_speed # Time it takes for the plates to move for every new layers
insertion_seed = 18

### NUMBER OF CELL PER DIRECTION
#... stuff need to be made
#...
### Next step

# Starts at 4, takes into account the reservoir, the build plate and the two seperators
solid_obj_id = 4
number_of_layers = number_of_layers + 2   # We add 2 for the layer -1 and 0
total_solid_objects = 4 + number_of_layers

# Finding the departure time of all the coaters
every_starting_time = np.empty(number_of_layers)
every_starting_time[0] = first_starting_time

# Times it takes for a coater to cross the domain
time_per_layer = lenght_domain / coater_speed

for i in range(1,number_of_layers):
    every_starting_time[i] = every_starting_time[i-1] + time_per_layer * 0.70
    # 0.7 is hard coded. This can be adjusted for simulation time optimization.
    # Normally, it is aroung 0.7 to 0.8 . The longer, the "saver" is it.

# Initializing the string
# First solid object is just for the leveling of the reservoir
end_time = (every_starting_time[solid_obj_id - 4] + time_per_layer + 0.05)
coater_function = f"if(t>= {every_starting_time[solid_obj_id - 4]:.5}, if(t<= {(every_starting_time[solid_obj_id - 4] + time_per_layer + 0.05):.5}, {coater_speed:.5}, 0), 0)"
all_coater= (f"   subsection solid object {solid_obj_id} \n"
             f"      subsection mesh \n"
             f"         set type                = gmsh \n"
             f"         set file name           = ./gmsh/Blade.msh\n"
             f"         set simplex             = true\n"
             f"         set initial refinement  = 0\n"
             f"         set initial translation = -0.001, 0, 0\n"
             f"      end\n"
             f"      subsection translational velocity\n"
             f"         set Function expression = {coater_function}; 0; 0\n"
             f"      end\n"
             f"   end\n")
solid_obj_id += 1

# This while loop is use to generate the translational velocity functions and to create
# the string that's gonna be replacing the "Blades" symbol in the parameter file.
while solid_obj_id - 4 < number_of_layers:
    end_time = (every_starting_time[solid_obj_id - 4] + time_per_layer + 0.05)
    coater_function = f"if(t>= {every_starting_time[solid_obj_id - 4]:.5}, if(t<= {(every_starting_time[solid_obj_id - 4] + time_per_layer + 0.05):.5}, {coater_speed:.5}, 0), 0)"

    coater_param = '''  subsection solid object ''' + str(solid_obj_id) + ''' 
    subsection mesh
      set type                = gmsh
      set file name           = ./gmsh/Blade.msh
      set simplex             = true
      set initial refinement  = 0
      set initial translation = -0.001, ''' + f"{delta_b_p:.5}" + ''', 0
    end
    subsection translational velocity
      set Function expression = ''' + coater_function +'''; 0; 0
    end
  end
   
'''
    solid_obj_id += 1
    all_coater = all_coater + coater_param


# Here, we create the displacement function for the reservoir plate and build plate.
reservoir_func   = str()
build_plate_func = str()

# The first one is special, since we want more powder on the first layer (2.2 mm). --> Layer number = -1
# It's the same as what you can find in the while loop right after.
t1 = every_starting_time[1] - plate_displacement_time * 4.  # We dont want big acceleration, thus de do the displacement over a longer period (4 - 1 = 3 times longer then a normal displacement)
t2 = every_starting_time[1] - plate_displacement_time * 1.
t3 = every_starting_time[1] + time_per_layer * 0.4 - plate_displacement_time * 2.
t4 = every_starting_time[1] + time_per_layer * 0.4 - plate_displacement_time * 1.

const_first_layer_extrusion  = first_layer_extrusion / ( (1/3) * (t2**3 - t1**3) + 0.5 * (t1+t2) * (t1**2 - t2**2) + (t1*t2) * (t2 - t1) )
const_first_layer_build_plate = delta_o              / ( (1/3) * (t4**3 - t3**3) + 0.5 * (t3+t4) * (t3**2 - t4**2) + (t3*t4) * (t4 - t3) )

reservoir_func   = reservoir_func   + f"if(t>= {t1:.5}, if(t<= {t2:.5}, { const_first_layer_extrusion:.5}  *(t - {t1:.5})*(t - {t2:.5}), "
build_plate_func = build_plate_func + f"if(t>= {t3:.5}, if(t<= {t4:.5}, {-const_first_layer_build_plate:.5}*(t - {t3:.5})*(t - {t4:.5}), "

for i in range(2,number_of_layers):
    t1 = every_starting_time[i] - plate_displacement_time * 2.                        # The reservoir moves...
    t2 = every_starting_time[i] - plate_displacement_time * 1.                        # ... and stops before the coater starts moving.
    t3 = every_starting_time[i] + time_per_layer * 0.4 - plate_displacement_time * 2. # The build-plate moves when the first coater is pass it...
    t4 = every_starting_time[i] + time_per_layer * 0.4 - plate_displacement_time * 1. # ... and stops before the second coater arrives.

    # This constant can be find by integrating the speed function we define. This speed function is a second degree polynomial.
    # The integral must be equal to the layer gap for the build plate.
    const = delta_n/ ( (1/3) * (t2**3 - t1**3) + 0.5 * (t1+t2) * (t1**2 - t2**2) + (t1*t2) * (t2 - t1) )

    reservoir_func   = reservoir_func   + f"if(t>= {t1:.5}, if(t<= {t2:.5}, {(const * other_layer_extrusion/delta_n):.5}*(t - {t1:.5})*(t - {t2:.5}), "
    build_plate_func = build_plate_func + f"if(t>= {t3:.5}, if(t<= {t4:.5}, {-const:.5}*(t - {t3:.5})*(t - {t4:.5}), "

# Ending the functions
reservoir_func   = reservoir_func   + "0) "
build_plate_func = build_plate_func + "0) "

for i in range(2, number_of_layers):
    reservoir_func   = reservoir_func   + ", 0) )"
    build_plate_func = build_plate_func + ", 0) )"

reservoir_func   = reservoir_func   + ", 0)"
build_plate_func = build_plate_func + ", 0)"

# Initial translation
initial_trans = f"{-(number_of_layers)        * 600E-6 - 2200E-6 :.5}"
y_min         = f"{-(number_of_layers + 0.5) * 600E-6 - 2200E-6 :.5}"


## File names ##
# Name of the .prm file created
CASE_PREFIX = f"lpbf_{delta_n*10**6:04.0f}_{coater_speed*10**3:04.0f}_{number_of_layers:02.0f}_{delta_o * 1e6 :03.0f}_{alloy}"

# Name of the output folder
output_folder = f"out_{delta_n*10**6:04.0f}_{coater_speed*10**3:04.0f}_{number_of_layers:02.0f}_{delta_o * 1e6 :03.0f}_{alloy}"

# Name of the restart files
restart_file = f"restart_{delta_n*10**6:04.0f}_{coater_speed*10**3:04.0f}_{number_of_layers:02.0f}_{delta_o * 1e6 :03.0f}_{alloy}"

# Name of the original .prm file being modified
PRM_FILE = 'lpbf_coating.prm'

# jinja2
PATH = os.getcwd()
templateLoader = jinja2.FileSystemLoader(searchpath=PATH)
templateEnv = jinja2.Environment(loader=templateLoader)
template = templateEnv.get_template(PRM_FILE)

# Replacing the symbols in the parameter file with the right expression
output_text = template.render(Coaters = all_coater,Total_solid_objects = str(total_solid_objects),
                              Reservoir_function = reservoir_func,
                              Build_plate_function = build_plate_func,
                              Initial_translation = initial_trans,
                              Y_min = y_min,
                              X_max = lenght_domain,
                              Out = output_folder,
                              End_time = end_time,
                              Seed = insertion_seed,
                              Diameters_List = diameters_list,
                              Restart = restart_file)

prm_file_name = CASE_PREFIX + ".prm"
output_file_path = os.path.join("./", prm_file_name)
with open(output_file_path, 'w') as f:
    f.write(output_text)
