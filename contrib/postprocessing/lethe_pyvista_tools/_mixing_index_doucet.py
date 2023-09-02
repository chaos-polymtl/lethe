import numpy as np
from scipy.linalg import eigh
from multiprocessing import Manager


# Calculate mixing index using the method by Doucet et al. (2008)
def mixing_index_doucet(self, reference_time_step = 0, use_cyl = False, increasing_index = False, normalize = True):
    """
    Calculates mixing index per time-step using the method by
    Doucet et al. (2008).
    J. Doucet, F. Bertrand, J. Chaouki. "A measure of mixing from
    Lagrangian tracking and its application to granular and fluid flow
    systems." Chemical Engineering Research and Design 86.12 (2008):
    1313-1321.

    Parameters:

    :param reference_time_step = 0     -> Time-step used as reference to
    calculate mixing index.

    :param use_cyl = False             -> Choose whether to use cylindrical or
    cartesian coordinates. If use_cyl = True, .point_cyl will be used
    (check get_cylindrical_coords method). Otherwise, cartesian .points are
    used.

    :param increasing_index = False    -> Choose whether the mixing index is
    increasing or decreasing with mixing. Doucet et al. (2008) uses a
    decreasing mixing index, however, most mixing indices increase with
    mixing.

    :param normalize = False           -> Choose whether the mixing index is
    normalized according to the mixing index of the reference_time_step.

    This method assigns the following attributes to the object:

    self.mixing_index           -> Normalized Doucet mixing index per
    time-step. The normalization is done using the mixing index at
    reference_time_step.

    self.mixing_eigenvector     -> Eigenvector associated to the
    mixing index.
    """

    # Apply method by  J. Doucet, F. Bertrand, J. Chaouki.
    # "A measure of mixing from Lagrangian tracking and its application to 
    # granular and fluid flow systems." Chemical Engineering Research and 
    # Design 86.12 (2008): 1313-1321.

    # Get cylindrical coordinates if requested and not previously
    # calculated

    if self.df_available:
        df = self.df[reference_time_step]
    else:
        df = self.get_df(reference_time_step)

    if use_cyl:
        self.get_cylindrical_coords()

        if not self.df_available:
            df = self.get_df(reference_time_step)

    # If cylindrical coordinates requested, assign points_cyl to reference
    # position, otherwise use cartesian
    if use_cyl:
        reference_position = df['points_cyl']
    else:
        reference_position = df.points

    # Get position of particles corresponding IDs
    id_keys = df["ID"]

    # Create list of mixing indices per time-step and array of eigenvectors
    self.mixing_index = Manager().list(np.empty(len(self.list_vtu)))
    self.mixing_eigenvector = Manager().list(np.empty(len(self.list_vtu)))

    # Loop through dataframes and find its mixing index
    global mixing_index_doucet_loop
    def mixing_index_doucet_loop(i):
        if self.df_available:
            df = self.df[i]
        else:
            df = self.get_df(i)

        # If cylindrical coordinates requested, assign points_cyl to current
        # position, otherwise use cartesian
        if use_cyl:
            i_position = df['points_cyl']
        else:
            i_position = df.points

        # Find indices of particles in different time-steps
        _, indices_i, indices_ref = np.intersect1d(df["ID"], id_keys, assume_unique = True, return_indices = True)

        # Calculate correlation matrix
        correlation_matrix = np.corrcoef(i_position[indices_i], reference_position[indices_ref], rowvar=False)[3:, :3]

        # Transpose matrix
        correlation_matrix_transpose = correlation_matrix.T

        # Multiply correlation and transposed correlation matrices
        M = np.matmul(correlation_matrix, correlation_matrix_transpose)

        # Find maximum eigenvalue and the eigenvector associated to it
        max_eigenvalue, assoc_eigenvectors = eigh(M, subset_by_index=[2, 2])
        max_eigenvalue = max_eigenvalue[0]

        # Store mixing index and associated eigenvector
        self.mixing_index[i] = max_eigenvalue
        self.mixing_eigenvector[i] = assoc_eigenvectors.flatten()

    self.parallel_run(mixing_index_doucet_loop, range(len(self.list_vtu)), tqdm_desc = "Calculating mixing index")

    # Normalize index
    if normalize:
        max_eigenvalue_reference = self.mixing_index[reference_time_step]
        self.mixing_index = np.divide(self.mixing_index, max_eigenvalue_reference)

    # Use increasing instead of decreasing index
    if increasing_index:
        self.mixing_index = 1 - self.mixing_index


