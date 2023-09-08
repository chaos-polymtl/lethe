import numpy as np


# Get cylindrical coordinates of each point of all dataframes
def get_cylindrical_coords(self, radial_components = "yz"):
    """
    Get cylindrical coordinates of points in self.df datasets

    Parameter:
    :param radial_components = "yz"           -> Cartesian directions of radial
    component.

    This method assigns the following attribute to the object:
    
    :return self.df[$TIME-STEP]['points_cyl'] -> Returns a .points like array with all
    points in cylindrical [radius, theta, height].
    """
    
    if self.has_cylindrical_coords:
        return

    # List of indices of radial components
    radial_indices = []

    # Add indices according to parameter radial_components
    if "x" in radial_components:
        radial_indices.append(0)
    if "y" in radial_components:
        radial_indices.append(1)
    if "z" in radial_components:
        radial_indices.append(2)

    # Kill process if radial_components have more or less than 2 coords
    if len(radial_components) != 2:
        print(f"radial_components has {len(radial_components)} axis")
        exit()
    
    # Find index other than the radial components
    z_index = [x for x in [0, 1, 2] if x not in radial_indices]

    # Loop through data
    global get_cylindrical_coords_loop
    def get_cylindrical_coords_loop(i):

        if self.df_available:
            df = self.df[i]
        else:
            df = self.get_df(i)

        # Get cartesian position
        cartesian = df.points

        # Calculate radial coord
        radius = np.sqrt(cartesian[:, radial_indices[0]]**2 + cartesian[:, radial_indices[1]]**2)

        # Calculate theta
        theta = np.arctan2(cartesian[:, radial_indices[1]], cartesian[:, radial_indices[0]])

        # Get z
        z = cartesian[:, z_index].flatten()

        # Store coordinates into points_cyl (same shape as .points)
        if self.df_available:
            self.df[i]['points_cyl'] = np.empty(df.points.shape)
            self.df[i]['points_cyl'][:, 0] = radius.tolist()
            self.df[i]['points_cyl'][:, 1] = theta
            self.df[i]['points_cyl'][:, 2] = z
        else:
            df['points_cyl'] = np.empty(df.points.shape)
            df['points_cyl'][:, 0] = radius.tolist()
            df['points_cyl'][:, 1] = theta
            df['points_cyl'][:, 2] = z
            df.save(f'{self.path_output}/{self.list_vtu[i]}')

    self.parallel_run(get_cylindrical_coords_loop, range(len(self.list_vtu)), tqdm_desc = "Getting cylindrical coords")

    self.has_cylindrical_coords = True