#############################################################################
"""Class of methods to post-process lethe results with pyvista"""
#############################################################################

# Import modules
import pyvista as pv
from tqdm import tqdm

# Define class:
class lethe_pyvista_tools():

    def __init__(self, case_path = '.', prm_file_name = '', pvd_name = '', prefix = '', first = 0, last = None,
                    step = 1, read_to_df = False, ignore_data = [], n_procs = None):
        """
        Constructor of post-processing object.
        
        The data follows the data models provided by PyVista.
        For further information, consult:
        https://docs.pyvista.org/user-guide/data_model.html
        
        The constructor parameters are (all accessible by self.$param):

        :param case_path = '.'         -> Path to the case, that is, the folder with
        the prm file. By default, the present folder.
        
        :param prm_file_name = ''      -> Name of the .prm file (including '.prm')
        with simulation setup.
        
        :param case_path               -> Path to the case.
        
        :param prm_file_name           -> Name of the prm_file.

        :param pvd_name                -> Name of the .pvd file containing the
        reference to Lethe data.

        :param prefix                  -> Prefix of the modified pvtu and pvd files. By default, 'mod_'.
        IMPORTANT!!!! If this parameter is empty, that is,  '',
        data will be written over the original pvtu and pvd files.

        :param first = 0               -> First time-step to be read into PyVista
        dataset.

        :param last  = None            -> Last time-step to be read into PyVista
        dataset.

        :param step  = 1               -> Step between datasets.

        :param read_to_df = False      -> Choose whether dataframes will be stored on
        RAM, that is, will be available on self.df list.

        :param ignore_data             -> List of data to be ignored when reading

        :param n_procs                 -> Number of processors used to run each function using parallel_run.
        If None, use the number of processors available.


        This method assigns the following attributes to the object:

        self.path_output        -> Returns the path to the output folder.

        self.prm_dict.get($PARAMETER)-> Returns the parameter value given a
        string with the parameter's name.
        
        self.pvd_name           -> Returns the name of the .pvd file.
        
        self.time_list          -> Returns the list of times corresponding to 
        datasets.

        self.list_vtu           -> Returns the list of names of .pvtu files.

        self.padding            -> Returns the padding of pvtu file numbering.

        """

        self.path_case = case_path
        self.prm_file = prm_file_name
        self.pvd_name = pvd_name
        self.ignore_data = ignore_data
        self.sorted = False
        self.has_neighbors = False
        self.has_cylindrical_coords = False
        self.padding = '0'

        if n_procs is None:
            from os import cpu_count
            self.n_procs = cpu_count()
        else:
            self.n_procs = n_procs

        if ".prm" not in self.prm_file:
            self.prm_file = self.prm_file + ".prm"
        
        # Read .prm file to dictionary
        # Create dictionary
        self.prm_dict = {}

        # Use .prm path as argument
        with open(self.path_case + '/' + self.prm_file) as file:

            # Loop through lines in .prm
            for line in file:

                # If the line has 'subsection' in it (and it is not commented)
                if 'subsection' in line and not '#' in line:

                    # Remove "subsection"
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

                        # If attribute already exists, create a list
                        # Otherwise, create key-value
                        if subsection_clean_line in self.prm_dict.keys():
                            if type(self.prm_dict[subsection_clean_line]) is list:
                                self.prm_dict[subsection_clean_line].append(clean_line[1])
                            
                            else:
                                self.prm_dict[subsection_clean_line] = [self.prm_dict[subsection_clean_line]]
                                self.prm_dict[subsection_clean_line].append(clean_line[1])
                    
                        else:
                            self.prm_dict[subsection_clean_line] = clean_line[1]
                    
                    else:
                    
                        # If attribute already exists, create a list
                        # Otherwise, create key-value
                        if clean_line[0] in self.prm_dict.keys():
                            if type(self.prm_dict[clean_line[0]]) is list:
                                self.prm_dict[clean_line[0]].append(clean_line[1])
                            
                            else:
                                self.prm_dict[clean_line[0]] = [self.prm_dict[clean_line[0]]]
                                self.prm_dict[clean_line[0]].append(clean_line[1])
                    
                        else:
                            self.prm_dict[clean_line[0]] = clean_line[1]

            print(f'Successfully constructed. To see the .prm dictionary, print($NAME.prm_dict)')
        
        # Define path where pvtu files are
        self.path_output = self.path_case + self.prm_dict['output path'].replace('.', '')
        
        # Read name of files in .pvd file        
        self.reader = pv.get_reader(f"{self.path_output}/{pvd_name}") 

        # Create list of pvd datasets
        pvd_datasets = self.reader.datasets

        # Create a list of time-steps
        self.time_list = self.reader.time_values

        # Create a list of all files' names
        list_vtu = [pvd_datasets[x].path for x in range(len(pvd_datasets))]

        # Remove duplicates
        list_vtu = list(dict.fromkeys(list_vtu))

        # Select data
        if last is None:
            self.list_vtu = list_vtu[first::step]
            self.time_list = self.time_list[first::step]
            self.first = first
            self.step = step
            self.last = len(self.time_list) - 1
            self.pvd_datasets = self.reader.datasets[first::step]
        else:
            self.list_vtu = list_vtu[first:last:step]
            self.time_list = self.time_list[first:last:step]
            self.first = first
            self.step = step
            self.last = last
            self.pvd_datasets = self.reader.datasets[first:last:step]

        if len(prefix) > 0:
            self.create_copy(prefix = prefix)

        # Boolean indicating that the dataframes are not stored in the
        # self.df object. If read_to_df = True is called, all data 
        # will be stored in self.df, thus, consuming a lot of RAM. 
        # Reading data into df can make the post-processing steps faster, since 
        # each step will be already available in self.df. However, this 
        # consumes a lot of RAM and for large simulations the tool will crash.
        # Alternatively, if read_to_df = False, all functions
        # will loop through the pvtu files and flush.
        self.df_available = False

        if read_to_df:
            # Create empty array to store results
            self.df = []

            # Read PVTU data
            n_pvtu = len(self.list_vtu)
            pbar = tqdm(total = n_pvtu, desc="Reading PVTU files")
            for i in range(len(self.list_vtu)):
                
                # Read dataframes from VTU files into df
                self.df.append(self.get_df)
                pbar.update(1)

            self.df_available = True

            print(f'Written .df[timestep] from timestep = 0 to timestep = {len(self.list_vtu)-1}')

    # IMPORT FUNCTIONS:

    # Apply multiprocessing imap
    from .parallel_run import parallel_run

    # Use prefix argument to create new pvtu and pvd
    from ._create_copy import create_copy

    # Return single pyvista dataset from list
    from ._get_df import get_df

    # Write modifications on each df to PVTU files
    from ._write_df_to_pvtu import write_df_to_pvtu

    # Sort all data given reference array
    from ._sort_by_array import sort_by_array

    # Create or modifie array
    from ._modify_array import modify_array

    # Get cylindrical coordinates of each point of all dataframes
    from ._get_cylindrical_coords import get_cylindrical_coords

    # Get neighbors of points
    from ._get_nearest_neighbors import get_nearest_neighbors

    # Calculate mixing index using Nearest Neighbors Method
    from ._mixing_index_nearest_neighbors import mixing_index_nearest_neighbors

    # Calculate mixing index using the method by Doucet et al. (2008)
    from ._mixing_index_doucet import mixing_index_doucet
