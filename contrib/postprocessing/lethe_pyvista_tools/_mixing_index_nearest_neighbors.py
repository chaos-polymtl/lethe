import numpy as np
from multiprocessing import Manager

# Calculate mixing index using Nearest Neighbors Method
def mixing_index_nearest_neighbors(self, n_neighbors = 15, reference_array = "particle_color", mixing_index_array_name = "mixing_index"):
    """
    Calculates mixing index per time-step using the Nearest
    Neighbors Method (NNM) by Godlieb et al. (2007).
    # Godlieb, W., N. G. Deen, and J. A. M. Kuipers. "Characterizing solids
    mixing in DEM simulations." cell 1 (2007).

    Parameters:
    :param n_neighbors = 15                        -> Number of neighbors to
    account for in the calculation of the mixing index.

    :param reference_array = "particle_color"      -> Name of the array containing
    the particle's type, that is, which group the particle is part of. For
    a better understanding, check the documentation of the modify_array
    method.

    :param mixing_index_array_name = "mixing_index"-> Name of the array assigned
    to each self.df. This array can be used to see the mixing index per
    point in the dataset (particle).

    This method assigns the following attributes to the object:

    self.mixing_index       -> Average mixing index per time-step.

    self.mixing_index_std   -> Standard deviation of the mixing index per
    time-step.

    self.df[$TIME-STEP]     -> Assign array to dataset named according to
    mixing_index_array_name. This array can be used in visual postprocessing
    softwares, such as ParaView. Check the write_vtu method of this module.
    """

    # Apply NNM by Godlieb et al. (2007)
    # Godlieb, W., N. G. Deen, and J. A. M. Kuipers.
    # "Characterizing solids mixing in DEM simulations." cell 1 (2007).

    # If neighbors is not an attribute of the dataframe

    self.n_neighbors = n_neighbors
    self.reference_array = reference_array
    self.mixing_index_array_name = mixing_index_array_name

    if self.df_available:
        df = self.df[0]
    else:
        df = self.get_df(0)

    self.get_nearest_neighbors(n_neighbors = n_neighbors)

    # Create empty list to store mixing_index per time-step
    self.mixing_index = Manager().list(np.empty(len(self.list_vtu)))
    self.mixing_index_std = Manager().list(np.empty(len(self.list_vtu)))

    # Loop through dataframes and find its mixing index
    global mixing_index_nearest_neighbors_loop
    def mixing_index_nearest_neighbors_loop(i):
        if self.df_available:
            df = self.df[i]
        else:
            df = self.get_df(i)

        # Find particles with different values for the reference array per
        # particle
        list_neighbor_reference_array = df[self.reference_array][df['neighbors']]
        n_equal_neighbors_per_particle = np.sum(np.equal(df[self.reference_array][:, None],
                                                         list_neighbor_reference_array), axis = 1)

        # Calculate mixing index per particle
        mixing_index_per_particle = 2*(1-(1/self.n_neighbors) * n_equal_neighbors_per_particle)

        # Create array of mixing index per particle
        if self.df_available:
            self.df[i][self.mixing_index_array_name] = mixing_index_per_particle
        else:
            df[self.mixing_index_array_name] = mixing_index_per_particle
            df.save(f'{self.path_output}/{self.list_vtu[i]}')

        mixing_index = np.mean(mixing_index_per_particle)
        mixing_index_std = np.std(mixing_index_per_particle)

        # Store mixing index
        self.mixing_index[i] = mixing_index
        self.mixing_index_std[i] = mixing_index_std

    self.parallel_run(mixing_index_nearest_neighbors_loop, range(len(self.list_vtu)), tqdm_desc = "Calculating mixing index")

