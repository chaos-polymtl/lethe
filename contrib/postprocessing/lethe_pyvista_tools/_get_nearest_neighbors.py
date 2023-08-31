from tqdm import tqdm
from sklearn.neighbors import KDTree

# Get neighbors of points
def get_nearest_neighbors(self, return_id = True, n_neighbors = 15):
    '''
    Get indices, distances, and "ID" (if requested) of nearest neighbors of 
    each point in self.df.

    Parameters:
    return_id = False         -> Decide whether ID is returned or not. If 
    True, but self.df does not have "ID", attribute neighbors_id is not
    assigned.

    This method assigns the following attributes to the object:
    
    self.df[$TIME-STEP].neighbors       -> Returns a lists with 
    indices of neighbors per point in dataset.

    self.df[$TIME-STEP].neighbors_dist  -> Returns a list with distances
    between neighbor points per point in dataset.

    self.df[$TIME-STEP].neighbors_id    -> Returns a list with "ID"
    of neighbor points per point in dataset.

    !!!IMPORTANT!!!
    
    This method uses KDTree to find neighbors.
    Details:
    https://scikit-learn.org/stable/modules/generated/sklearn.neighbors.KDTree.html

    If library sklearn is missing, the method will not work.
    To install sklearn, run the following line in your terminal:
    $ pip install scikit-learn
    or
    $ pip3 install scikit-learn
    At the moment of the implementation, sklearn version was: 1.2.1
    '''

    # Loop through dataframes to search for neighbors
    pbar = tqdm(total = len(self.list_vtu), desc = "Finding neighbors")
    for i in range(len(self.list_vtu)):

        if self.df_available:
            df = self.df[i]
        else:
            df = self.get_df(i)

        # Create a tree from points
        tree = KDTree(df.points)

        # Get the distance and the indices of the n_neighbors neighbors
        # It is important to note that the closest neighbor is going to
        # be the point itself, so we ask for n_neighbors + 1
        dist, indices = tree.query(df.points, k = n_neighbors+1)

        # Remove itself from indices and dist for all points
        indices = indices[:, 1:]
        dist = dist[:, 1:]

        # Add neighbors_id, neighbors indices, and neighbors distances
        # to each dataframe
        if self.df_available:
            if return_id and hasattr(df, "ID"):
                self.df[i]["neighbors_id"] = self.df[i]["ID"][indices]
            self.df[i]["neighbors"] = indices
            self.df[i]["neighbors_dist"] = dist
        else:
            if return_id and hasattr(df, "ID"):
                df["neighbors_id"] = df["ID"][indices]
            df["neighbors"] = indices
            df["neighbors_dist"] = dist
            df.save(f'{self.path_output}/{self.list_vtu[i]}')

        pbar.update(1)