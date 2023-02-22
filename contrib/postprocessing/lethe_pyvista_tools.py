#############################################################################
"""Class of methods to post-process lethe results with pyvista"""
#############################################################################

# Import modules
import os
import numpy as np
import pyvista as pv
from tqdm import tqdm
from operator import itemgetter

# Define class:
class lethe_pyvista_tools():

    def __init__(self, case_path, prm_file_name):

        self.path_case = case_path
        self.prm_file = prm_file_name

        # Read .prm file to dictionary
        # Create dictionary
        self.prm_dict = {}

        # Use .prm path as argument
        with open(self.path_case + '/' + self.prm_file) as file:
            # Loop through lines in .prm
            for line in file:
                # If the line has 'subsection' in it (and it is not commented)
                if 'subsection' in line and not '#' in line:
                    # Remove "subsetction"
                    subsection_clean_line = line.replace('subsection', '')
                    # Clean line from spaces and assign key-value
                    subsection_clean_line = subsection_clean_line.strip()
                # Else, if the line has 'set' in it (and it is not commented)
                elif 'set' in line and not '#' in line:
                    # Remove "set" from string "line"
                    clean_line = line.replace('set', '')
                    # Split the string in [variable, value]
                    clean_line = clean_line.split('=')
                    # Clean line from spaces
                    for element in range(len(clean_line)):
                        clean_line[element] = clean_line[element].strip()
                    # Convert values to float when possible
                    try:
                        clean_line[1] = float(clean_line[1])
                    except:
                        pass
                    # Define [variable, value] as key and value in the
                    # dictionary
                    # If 'set' is a 'Function expression' or 'type'
                    if clean_line[0] == 'Function expression' or clean_line[0] == 'type':
                        # Attribute the name of the subsection above
                        self.prm_dict[subsection_clean_line] = clean_line[1]
                    else:
                        # Otherwise, attribute the set name
                        self.prm_dict[clean_line[0]] = clean_line[1]

            print(f'Successfully constructed. To see the .prm dictionary, print($NAME.prm_dict)')
        # Define path where vtu files are
        self.path_output = self.path_case + self.prm_dict['output path'].replace('.', '')

    # Read fluid or particle information from vtu files
    def read_lethe_to_pyvista(self, pvd_name, first = 0, last = None, interval = 1):

        self.pvd_name = pvd_name
        # Read name of files in .pvd file        
        reader = pv.get_reader(f"{self.path_output}/{pvd_name}") 

        # Create list of pvd datasets
        self.pvd_datasets = reader.datasets

        # Create a list of time-steps
        self.time_list = reader.time_values

        # Create a list of all files' names
        self.list_vtu = [self.pvd_datasets[x].path for x in range(len(self.pvd_datasets))]
        self.list_vtu = [x.replace(".pvtu", ".0000.vtu") for x in self.list_vtu]

        if last == None:
            self.list_vtu = self.list_vtu[first::interval]
            self.time_list = self.time_list[first::interval]
            self.first = first
            self.interval = interval
            self.last = len(self.time_list) - 1
        else:
            self.list_vtu = self.list_vtu[first:last:interval]
            self.time_list = self.time_list[first:last:interval]
            self.first = first
            self.interval = interval
            self.last = last

        # Create empty list to store results
        self.df = []

        # Read VTU data
        N_vtu = len(self.list_vtu)
        pbar = tqdm(total = N_vtu, desc="Reading VTU files")
        for i in range(len(self.list_vtu)):
            # Read dataframes from VTU files into df
            self.df.append(pv.read(f"{self.path_output}/{self.list_vtu[i]}"))
            pbar.update(1)

        print(f'Written .df[timestep] from timestep = 0 to timestep = {len(self.list_vtu)-1}')

    # Write modifications on each df to VTU files
    def write_vtu(self, prefix = "mod_"):

        # List of paths among read data
        read_files_path_list = [self.pvd_datasets[x].path for x in range(len(self.pvd_datasets))]

        # Write modified PVD to match new VTU files
        with open(f'{self.path_output}/{self.pvd_name}') as pvd_in:
            with open(f'{self.path_output}/mod_{self.pvd_name}', 'w') as pvd_out:
                for line in pvd_in:
                    
                    # If line refers to a dataset
                    if "vtu" in line:

                        # For all read files
                        for path in read_files_path_list:

                            # If line matches one of the files
                            if path in line:
                                line = line.replace('.pvtu', '.0000.vtu')
                                line = line.replace('file="', f'file="{prefix}')
                                pvd_out.write(line)
                                read_files_path_list.remove(path)
                                pass
                    
                    # Write config lines
                    else:
                        pvd_out.write(line)
        
        # Write modified VTU file
        N_vtu = len(self.df)
        pbar = tqdm(total = N_vtu, desc="Writting new VTU and PVD files")
        for i in range(len(self.df)):
            self.df[i].save(f'{self.path_output}/{prefix}{self.list_vtu[i]}')
            pbar.update(1)


        print(f"Modified .vtu and .pvd files with prefix {prefix} successfully written")

    # Sort all data given reference array 
    def sort_by_array(self, reference_array_name):
        
        pbar = tqdm(total = len(self.time_list), desc = f"Sorting dataframe by {reference_array_name}")
        for i in range(len(self.time_list)):
            self.df[i].points = self.df[i].points[self.df[i][reference_array_name].argsort()]
            for name in self.df[0].array_names:
                self.df[i][name] = self.df[i][name][self.df[i][reference_array_name].argsort()]
            pbar.update(1)

    # Creates or modifies array
    def array_modifier(self, reference_array_name = "ID", array_name = "new_array", restart_array = False,  condition = "", array_values = 0, standard_value = 0, reference_time_step = 0, time_dependent = False):
        '''
        Parameters are:

        reference_array_name = "ID"        -> array to be used as reference to
        create or modify the other. All arrays will be sorted and written according 
        to this one.

        array_name           = "new_array" -> name of the new array. If there is an
        array with the same name, it will be rewritten according to the other
        arguments.

        restart_array        = False       -> if True, zeroes the entire array before
        modifying it. If you want to modify part of the array keeping the rest
        intact, set it to False

        condition            = ""          -> takes a string and uses it in an if
        condition to modify the array. Variables accepted include x, y, z, u, v, w,
        t, and any other array.It also accepts a combination of them, such as:
        "x*w**2 + t > 2"

        array_values         = 0           -> new values to the array. This argument
        accepts a single value (which will be repeated to all data respecting the
        given condition), an numpy array or python list (with the same len of all
        other arrays), or a string such as "2*x + t" (working just like the condition
        argument)

        standard_value       = 0           -> if restart array is True, the
        standard_value will be the one plugged to the entire array before modifying
        it.

        reference_time_step  = 0           -> reference time step to which the
        modification will be applied. The others will follow this given one.

        time_dependent       = False       -> the modifier can be time dependent or
        not. If set True, the condition will be tested to each of the time-steps,
        while if False, it will be applied using the reference_time_step instead, and
        the modification will be just replicated to the other time steps
        '''

        print("Generating array based on condition and array_value")

        # Sort all data by reference_array_name
        print(f"Sort array by {reference_array_name}")
        self.sort_by_array(reference_array_name)

        # Time array
        t = self.time_list

        # Create list of array names
        # This step is necessary to allow the usage of
        # the variables x, y, z, u, v, w, t, f_x, f_y, and f_z
        # in the condition argument
        array_names = self.df[0].array_names
        array_names.append("x")
        array_names.append("y")
        array_names.append("z")
        array_names.append("u")
        array_names.append("v")
        array_names.append("w")
        array_names.append("f_x")
        array_names.append("f_y")
        array_names.append("f_z")
        array_names.append("t")

        # Restart array if asked or if array does not exist
        # If the array is restarted or created, the standard_value will be
        # assigned to the entire array.
        # If restart_array is set to False and the array exists,
        # The previous values in it will be preserved.
        # This can be used to apply multiple conditions without affecting
        # Previous modifications, for example. 
        if restart_array == True or array_name not in array_names:
            # Create array if does not exist
            new_array = np.repeat(standard_value, len(self.df[reference_time_step][reference_array_name]))
            print(f"Creating array '{array_name}' with standard_value {standard_value}")

            # Push array to all pyvista arrays
            pbar = tqdm(total = len(self.df), desc = f"Creating array: {array_name}")
            for i in range(len(self.df)):
                self.df[i][array_name] = np.repeat(standard_value, len(self.df[i][reference_array_name]))
                pbar.update(1)

        else:
            # Reading array from reference timestep
            print("Reading previous array")
            new_array = self.df[reference_time_step][array_name]


        # Create a list of array names that are used either in
        # "conditions" or in "array_values"
        new_variables = set([])

        # Prepare "condition" and "array_value" for elementwise loop
        # Note that "k" is used here because it is the specific counter
        # that will be used for testing the "condition" further
        for name in array_names:
            if name in condition:
                condition = condition.replace(name, name + "[k]")

                # If one of the variables used in "condition"
                # is a pyvista array, create a list with the
                # name of the variable for further manipulation
                if name in self.df[0].array_names:
                    exec(f"global {name}; {name} = self.df[reference_time_step][name]")
                    new_variables.add(name)

        if type(array_values) == type(str()):
            for name in array_names:
                if name in array_values:
                    array_values = array_values.replace(name, name + "[k]")

                    # If one of the variable used in "array_value"
                    # is a pyvista array, create a list with the
                    # name of the variable for further manipulation
                    if name in self.df[0].array_names:
                        exec(f"global {name}; {name} = self.df[reference_time_step][name]")
                        new_variables.add(name)

        # If results vary with time,
        # the condition and array_values will be applied
        # to all time steps
        if time_dependent:
            ("Creating time-dependent array:")
            pbar = tqdm(total = len(self.time_list), desc = f"Looping through time-steps")
            for i in range(len(self.df)):
                # Assign velocities and positions to variables using the ith
                # time step
                exec(f"global x; x = self.df[i].points[:, 0]")
                exec(f'global y; y = self.df[i].points[:, 1]')
                exec(f'global z; z = self.df[i].points[:, 2]')


                # In case velocity is written with caps V or v
                if "velocity" in self.df[0].array_names:
                    exec(f'global u; u = self.df[i]["velocity"][:, 0]')
                    exec(f'global v; v = self.df[i]["velocity"][:, 1]')
                    exec(f'global w; w = self.df[i]["velocity"][:, 2]')


                elif "Velocity" in self.df[0].array_names:
                    exec(f'global u; u = self.df[i]["Velocity"][:, 0]')
                    exec(f'global v; v = self.df[i]["Velocity"][:, 1]')
                    exec(f'global w; w = self.df[i]["Velocity"][:, 2]')

                # In case of FemForce
                if "FemForce" in self.df[0].array_names:
                    exec(f'global f_x; f_x = self.df[i]["FemForce"][:, 0]')
                    exec(f'global f_y; f_y = self.df[i]["FemForce"][:, 1]')
                    exec(f'global f_z; f_z = self.df[i]["FemForce"][:, 2]')

                # Update lists used either in "condition" or "array_value":
                for variable in new_variables:
                    exec(f"{variable} = self.df[i][variable]")

                # Reading array from reference timestep
                new_array = self.df[i][array_name]

                # Fill new_array with array_value
                for k in range(len(new_array)):
                    if eval(condition):
                        if type (array_values) == type(int(1)):
                            new_array[k] = array_values
                        elif type(array_values) == type(np.array([])) or type(array_values) == type([]):
                            new_array[k] = array_values[k]
                        else:
                            new_array[k] = eval(array_values)

                # Assign new_array to pyvista dataframe
                self.df[i][array_name] = new_array
                pbar.update(1)

        # If not time dependent, the condition and array_values will be applied
        # at the reference_time_step.
        # This is very useful if you want to track a group of particles in a
        # DEM simulation. The value of 1 can be assigned to this group of
        # particles
        # according to a given condition (position, for example) while 0 can be
        # given to the others,
        # This way, the groups of particles will be colored in the
        # reference_time_step and will keep this color regardless of time.
        else:
            print(f"Creating array based on time-step number: {reference_time_step}")
            print(f"Corresponding time: {self.time_list[reference_time_step]}")
            # Assign velocities and positions to variables using
            # reference_time_step
            exec(f'global x; x = self.df[reference_time_step].points[:, 0]')
            exec(f'global y; y = self.df[reference_time_step].points[:, 1]')
            exec(f'global z; z = self.df[reference_time_step].points[:, 2]')


            # In case velocity is written with caps V or v
            if "velocity" in self.df[0].array_names:
                exec(f'global u; u = self.df[reference_time_step]["velocity"][:, 0]')
                exec(f'global v; v = self.df[reference_time_step]["velocity"][:, 1]')
                exec(f'global w; w = self.df[reference_time_step]["velocity"][:, 2]')


            elif "Velocity" in self.df[0].array_names:
                exec(f'global u; u = self.df[reference_time_step]["Velocity"][:, 0]')
                exec(f'global v; v = self.df[reference_time_step]["Velocity"][:, 1]')
                exec(f'global w; w = self.df[reference_time_step]["Velocity"][:, 2]')

            # In case of FemForce
            if "FemForce" in self.df[0].array_names:
                exec(f'global f_x; f_x = self.df[reference_time_step]["FemForce"][:, 0]')
                exec(f'global f_y; f_y = self.df[reference_time_step]["FemForce"][:, 1]')
                exec(f'global f_z; f_z = self.df[reference_time_step]["FemForce"][:, 2]')

            # Fill new_array with array_value
            print(f"Creating new array named: {array_name}")
            for k in range(len(new_array)):
                if eval(condition):
                    if type (array_values) == type(int(1)):
                        new_array[k] = array_values
                    elif type(array_values) == type(np.array([])) or type(array_values) == type([]):
                        new_array[k] = array_values[k]
                    else:
                        new_array[k] = eval(array_values)

            # Assign new_array to pyvista dataframe
            self.df[reference_time_step][array_name] = new_array

            # Create dictionary (map) based on reference_array
            reference_time_step_dict = dict(zip(self.df[reference_time_step][reference_array_name], self.df[reference_time_step][array_name]))
            
            key_list = self.df[reference_time_step][reference_array_name]

            # Use the same values for all time steps
            # Note that "reference_array_name" is used as criterium here
            # for sorting purposes, and that it can be changed
            # according to the user by changin the parameter
            # "reference_array_name" to any other array name in the original
            # pyvista arrays
            # (self.df[0].array_names, for example)
            pbar = tqdm(total = len(self.df), desc = f"Assigning {array_name} to dataframes")
            for i in range(len(self.df)):

                # Find elements in common in current and reference arrays
                keys, indices, _ = np.intersect1d(self.df[i][reference_array_name], key_list, assume_unique = True, return_indices = True)
                self.df[i][array_name][indices] = itemgetter(*keys)(reference_time_step_dict)

                pbar.update(1)
    
