# To use the lethe_pyvista_tools, you need to have python 3
# installed in your computer

# The modules necessary to run lethe pyvista tools are:
# Pandas: pip install pandas
# PyVista: pip install pyvista

# To use lethe_pyvista_tools, append the path to the /contrib/postprocessing
# folder in Lethe to sys 
import sys

# Path to the module
path_to_module = '$LETHE_PATH/contrib/postprocessing/'
sys.path.append(path_to_module)

# or put the "lethe_pyvista_tools.py" file inside
# the same directory as your python script and procceed as follows

# This line imports all lethe_pyvista_tools functionalities
from lethe_pyvista_tools import *

# This script prints out the content of your prm file as dictionary
# To run lethe_pyvista_tools you need to specify the path to your
# case and the name of the .prm file
example = lethe_pyvista_tools('PATH TO YOUR CASE', 'NAME_OF_YOUR_PARAMETERS_FILE.prm')

print('This is the dictionary of your .prm file:')
print(example.prm_dict)
print('To print out any value inside the dictionary, ask for it using ["parameter_name"] right after .prm_dict variable')

print('The path to the case can be seen using: example.case_path')

# To read the data to pyvista dataframe, use the following with
# the .pvd file as argument
example.read_lethe_to_pyvista('NAME_OF_YOUR_PVD_FILE.pvd')

# The read_lethe_to_pyvista method writes out the attributes .time_list,
# .list_vtu and reads the '.vtu' files inside the pointed folder as pyvista
# dataframes.
print('List of all .vtu: ')
print(example.list_vtu)
print('Time list, if transient: ')
print(example.time_list)

# Each .vtu file will correspond to a dataframe named df, such that the first
# vtu can be accessed through .df[0], the second .df[1], and so on.
print(example.df[0])

# This should print out the name of the arrays in the first vtu file of your
# case
print('Name of the arrays in your pyvista dataframe: ')
print(example.df[0].array_names)

# This function sorts all arrays according to a given one array name.
# For particle data, for example, the particles can be sorted according
# to their ID using:
example.sort_by_array("ID")

# It is also possible to modify or create new arrays using
example.array_modifier()
# The arguments (and their default values) are:

# reference_array_name = "ID"        -> array to be used as reference to
# create or modify the other. All arrays will be sorted and written according 
# to this one.

# array_name           = "new_array" -> name of the new array. If there is an
# array with the same name, it will be rewritten according to the other
# arguments.

# restart_array        = False       -> if True, zeroes the entire array before
# modifying it. If you want to modify part of the array keeping the rest
# intact, set it to False

# condition            = ""          -> takes a string and uses it in an if
# condition to modify the array. Variables accepted include x, y, z, u, v, w,
# t, and any other array.It also accepts a combination of them, such as:
# "x*w**2 + t > 2"

# array_values         = 0           -> new values to the array. This argument
# accepts a single value (which will be repeated to all data respecting the
# given condition), an numpy array or python list (with the same len of all
# other arrays), or a string such as "2*x + t" (working just like the condition
# argument)

# standard_value       = 0           -> if restart array is True, the
# standard_value will be the one plugged to the entire array before modifying
# it.

# reference_time_step  = 0           -> reference time step to which the
# modification will be applied. The others will follow this given one.

# time_dependent       = False       -> the modifier can be time dependent or
# not. If set True, the condition will be tested to each of the time-steps,
# while if False, it will be applied using the reference_time_step instead, and
# the modification will be just replicated to the other time steps


# Other modifications can be made without using the array_modified function, 
# using the same structure as in PyVista documentation
# (https://docs.pyvista.org/)

# To write those, the following can be used:
example.write_vtu()

# For further information about PyVista, refer to https://docs.pyvista.org/