import pyvista as pv


# Return single pyvista dataset from list
def get_df(self, time_step):
    """
    Reads and returns data from .pvtu file.

    :param time_step -> Time-step number.
    :return: PyVista array with data.
    """

    # Get reader for the PVTU file
    pvtu_reader = pv.get_reader(f"{self.path_output}/{self.list_vtu[time_step]}")

    # Ignore selected data in order to reduce RAM usage
    for data in self.ignore_data:
        pvtu_reader.disable_point_array(data)

    data = pvtu_reader.read()
    data = data.cast_to_unstructured_grid()

    return data