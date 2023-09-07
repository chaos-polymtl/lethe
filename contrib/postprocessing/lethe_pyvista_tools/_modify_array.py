import numpy as np
from operator import itemgetter


# Creates or modifies array
def modify_array(self, reference_array_name = "ID", array_name = "new_array", restart_array = False,  condition = "",
                 array_values = 0, standard_value = 0, reference_time_step = 0, time_dependent = False):
    """
    Creates or modifies array

    Parameters are:
    :param reference_array_name = "ID"        -> array to be used as reference to
    create or modify the other. All arrays will be sorted and written according to this one.

    :param array_name           = "new_array" -> name of the new array. If there is an
    array with the same name, it will be rewritten according to the other
    arguments.

    :param restart_array        = False       -> if True, zeroes the entire array before
    modifying it. If you want to modify part of the array keeping the rest
    intact, set it to False

    :param condition            = ""          -> takes a string and uses it in an if
    condition to modify the array. Variables accepted include x, y, z, u, v, w,
    t, and any other array.It also accepts a combination of them, such as:
    "x*w**2 + t > 2"

    :param array_values         = 0           -> new values to the array. This argument
    accepts a single value (which will be repeated to all data respecting the
    given condition), an numpy array or python list (with the same length as all

    other arrays), or a string such as "2*x + t" (working just like the condition
    argument)

    :param standard_value       = 0           -> if restart array is True, the
    standard_value will be the one plugged to the entire array before modifying
    it.

    :param reference_time_step  = 0           -> reference time step to which the
    modification will be applied. The others will follow this given one.

    :param time_dependent       = False       -> the modifier can be time dependent or
    not. If set True, the condition will be tested to each of the time-steps,
    while if False, it will be applied using the reference_time_step instead, and
    the modification will be just replicated to the other time steps
    """

    print("Generating array based on condition and array_value")

    # Sort all data by reference_array_name
    self.sort_by_array(reference_array_name)

    # Time array
    t = self.time_list

    # Create list of array names
    # This step is necessary to allow the usage of
    # the variables x, y, z, u, v, w, t, f_x, f_y, and f_z
    # in the condition argument
    if self.df_available:
        df = self.df[0]
    else:
        df = self.get_df(0)
    array_names = df.array_names
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
    if restart_array or array_name not in df.array_names:
        # Create array if it does not exist
        new_array = np.repeat(standard_value, len(df[reference_array_name]))
        print(f"Creating array '{array_name}' with standard_value {standard_value}")

        # Push array to all pyvista arrays
        global create_array
        def create_array(i):
            if self.df_available:
                self.df[i][array_name] = np.repeat(standard_value, len(self.df[i][reference_array_name]))
            
            else:
                df = self.get_df(i)
                df[array_name] = np.repeat(standard_value, len(df[reference_array_name]))
                df.save(f'{self.path_output}/{self.list_vtu[i]}')

        self.parallel_run(create_array, range(len(self.list_vtu)), tqdm_desc = f"Creating array {array_name}")

    else:
        if self.df_available:
            # Reading array from reference timestep
            print("Reading previous array")
            new_array = self.df[reference_time_step][array_name]

        else:
            print("Reading previous array")
            df_reference = self.get_df(reference_time_step)
            new_array = df_reference[array_name]

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
            if name in df.array_names:
                if self.df_available:
                    exec(f"global {name}; {name} = self.df[reference_time_step][name]")
                else:
                    df_reference = self.get_df(reference_time_step)
                    exec(f"global {name}; {name} = df_reference[name]")
                new_variables.add(name)

    if type(array_values) == type(str()):
        for name in array_names:
            if name in array_values:
                array_values = array_values.replace(name, name + "[k]")

                # If one of the variable used in "array_value"
                # is a pyvista array, create a list with the
                # name of the variable for further manipulation
                if name in df.array_names:
                    if self.df_available:
                        exec(f"global {name}; {name} = self.df[reference_time_step][name]")
                    else:
                        df_reference = self.get_df(reference_time_step)
                        exec(f"global {name}; {name} = df_reference[name]")
                    new_variables.add(name)

    # If results vary with time,
    # the condition and array_values will be applied
    # to all time steps
    global modify_array_loop
    if time_dependent:
        def modify_array_loop(i):
            # Assign velocities and positions to variables using the ith
            # time step
            if self.df_available:
                df = self.df[i]

            else:
                df = self.get_df(i)

            x = df.points[:, 0]
            y = df.points[:, 1]
            z = df.points[:, 2]


            # In case velocity is written with caps V or v
            if "velocity" in array_names:
                u = df["velocity"][:, 0]
                v = df["velocity"][:, 1]
                w = df["velocity"][:, 2]


            elif "Velocity" in array_names:
                u = df["Velocity"][:, 0]
                v = df["Velocity"][:, 1]
                w = df["Velocity"][:, 2]

            # In case of FemForce or fem_force
            if "FemForce" in array_names:
                f_x = df["FemForce"][:, 0]
                f_y = df["FemForce"][:, 1]
                f_z = df["FemForce"][:, 2]

            elif "fem_force" in array_names:
                f_x = df["fem_force"][:, 0]
                f_y = df["fem_force"][:, 1]
                f_z = df["fem_force"][:, 2]

            # In case of fem_torque
            if "fem_torque" in array_names:
                t_x = df["fem_torque"][:, 0]
                t_y = df["fem_torque"][:, 1]
                t_z = df["fem_torque"][:, 2]

            # Update lists used either in "condition" or "array_value":
            for variable in new_variables:
                exec(f"{variable} = df[variable]")

            # Reading array from reference timestep
            new_array = df[array_name]

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
            df[array_name] = new_array
            
            if self.df_available:
                self.df[i] = df
            else:
                df.save(f'{self.path_output}/{self.list_vtu[i]}')

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
        
        if self.df_available:
            df_reference = self.df[reference_time_step]

        else:
            df_reference = self.get_df(reference_time_step)

        x = df_reference.points[:, 0]
        y = df_reference.points[:, 1]
        z = df_reference.points[:, 2]


        # In case velocity is written with caps V or v
        if "velocity" in array_names:
            u = df_reference["velocity"][:, 0]
            v = df_reference["velocity"][:, 1]
            w = df_reference["velocity"][:, 2]


        elif "Velocity" in array_names:
            u = df_reference["Velocity"][:, 0]
            v = df_reference["Velocity"][:, 1]
            w = df_reference["Velocity"][:, 2]

        # In case of FemForce or fem_force
        if "FemForce" in array_names:
            f_x = df_reference["FemForce"][:, 0]
            f_y = df_reference["FemForce"][:, 1]
            f_z = df_reference["FemForce"][:, 2]

        elif "fem_force" in array_names:
            f_x = df_reference["fem_force"][:, 0]
            f_y = df_reference["fem_force"][:, 1]
            f_z = df_reference["fem_force"][:, 2]

        # In case of fem_torque
        if "fem_torque" in array_names:
            t_x = df_reference["fem_torque"][:, 0]
            t_y = df_reference["fem_torque"][:, 1]
            t_z = df_reference["fem_torque"][:, 2]

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
        if self.df_available:
            self.df[reference_time_step][array_name] = new_array
        else:
            df_reference[array_name] = new_array
            df_reference.save(f'{self.path_output}/{self.list_vtu[reference_time_step]}')

        # Create dictionary (map) based on reference_array
        reference_time_step_dict = dict(zip(df_reference[reference_array_name], df_reference[array_name]))
        
        key_list = df_reference[reference_array_name]

        # Use the same values for all time steps
        # Note that "reference_array_name" is used as criterium here
        # for sorting purposes, and that it can be changed
        # according to the user by changing the parameter
        # "reference_array_name" to any other array name in the original
        # pyvista arrays

        def modify_array_loop(i):
            # Find elements in common in current and reference arrays
            if self.df_available:
                df = self.df[i]
            else:
                df = self.get_df(i)
            
            keys, indices, _ = np.intersect1d(df[reference_array_name], key_list, assume_unique = True, return_indices = True)

            if self.df_available:
                self.df[i][array_name][indices] = itemgetter(*keys)(reference_time_step_dict)
            else:
                df[array_name][indices] = itemgetter(*keys)(reference_time_step_dict)
                df.save(f'{self.path_output}/{self.list_vtu[i]}')

    self.parallel_run(modify_array_loop, range(len(self.list_vtu)), tqdm_desc = "Assigning array")
