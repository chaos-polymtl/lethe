import pyvista as pv


# Return single pyvista dataset from list
def get_df(self, time_step):
    """
    Reads and returns data from .vtu file.

    :param time_step -> Time-step number.
    :return: PyVista array with data.
    """

    # Get reader for the VTU file
    vtu_reader = pv.get_reader(f"{self.path_output}/{self.list_vtu[time_step]}")

    # Ignore selected data in order to reduce RAM usage
    for data in self.ignore_data:
        vtu_reader.disable_point_array(data)

    return vtu_reader.read()