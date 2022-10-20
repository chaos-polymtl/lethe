#############################################################################
"""Class of methods to post-process lethe results with pyvista"""
#############################################################################

#Import modules
import numpy as np
import pandas as pd
import pyvista as pv
from tqdm import tqdm
from multiprocessing import Pool, cpu_count

#Define class:
class Lethe_pyvista_tools():

    def __init__(self, case_path, prm_file_name):
        
        self.path_case = case_path
        self.prm_file = prm_file_name
    
        #Read .prm file to dictionary
        #Create dictionary
        self.prm_dict = {}
        

        #Use .prm path as argument
        with open(self.path_case + '/' + self.prm_file) as file:
            #Loop trhough lines in .prm
            for line in file:
                #If the line has 'subsection' in it (and it is not commented)
                if 'subsection' in line and not '#' in line:
                    #Remove "subsetction"
                    subsection_clean_line = line.replace('subsection', '')
                    #Clean line from spaces and assign key-value
                    subsection_clean_line = subsection_clean_line.strip()
                #Else, if the line has 'set' in it (and it is not commented)
                elif 'set' in line and not '#' in line:
                    #Remove "set" from string "line"
                    clean_line = line.replace('set', '')
                    #Split the string in [variable, value]
                    clean_line = clean_line.split('=')
                    #Clean line from spaces
                    for element in range(len(clean_line)):
                        clean_line[element] = clean_line[element].strip()
                    #Convert values to float when possible
                    try:
                        clean_line[1] = float(clean_line[1])
                    except:
                        pass
                    #Define [variable, value] as key and value in the dictionary
                    #If 'set' is a 'Function expression' or 'type'
                    if clean_line[0] == 'Function expression' or clean_line[0] == 'type':
                        #Attribute the name of the subsection above
                        self.prm_dict[subsection_clean_line] = clean_line[1]
                    else:
                        #Otherwise, attribute the set name
                        self.prm_dict[clean_line[0]] = clean_line[1]
        
            print(f'Successfully constructed. To see the .prm dictionary, print($NAME.prm_dict)')
        #Define path where vtu files are
        self.path_output = self.path_case + self.prm_dict['output path'].replace('.', '')

    def reader(self, i = 0):
        exec(f'self.df_{i} = pv.read(f\'{self.path_output}/{self.list_vtu[i]}\')')
        #df = pv.read(f'{self.path_output}/{self.list_vtu[i]}')
        #return df

    #Read fluid information from vtu files
    def read_lethe_to_pyvista(self, pvd_name, first = 0, last = None, interval = 1):
        
        #Read name of files in .pvd file
        files = pd.read_csv(f'{self.path_output}{pvd_name}',sep='"',skiprows=6, usecols=[1, 5], names = ['time', 'vtu'])
        #clean data from NaN's
        files = files.dropna()
        #Create a list of time-steps
        self.time_list = files['time'].tolist()
        #Create a list of all files' names
        self.list_vtu = files['vtu'].tolist()
        #Format files' names
        self.list_vtu = [i.replace('.pvtu', '.0000.vtu') for i in self.list_vtu]

        if last == None:
            self.list_vtu = self.list_vtu[first::interval]
            self.time_list = self.time_list[first::interval]
        else:
            self.list_vtu = self.list_vtu[first:last:interval]
            self.time_list = self.time_list[first:last:interval]

        #Read VTU data
        N_vtu = len(self.list_vtu)
        pbar = tqdm(total = N_vtu, desc="Reading VTU files")
        for i in range(len(self.list_vtu)):
            #Read DF from VTU files
            self.reader(i)#exec(f'self.df_{i} = pv.read(f\'{self.path_output}/{self.list_vtu[i]}\')')
            pbar.update(1)
        
        print(f'Written .df_timestep from timestep = 0 to timestep = {len(self.list_vtu)-1}')
    
    def read_lethe_to_pyvista_parallel(self, pvd_name, first = 0, last = None, interval = 1, processors = cpu_count()):
        #Read name of files in .pvd file
        files = pd.read_csv(f'{self.path_output}{pvd_name}',sep='"',skiprows=6, usecols=[1, 5], names = ['time', 'vtu'])
        #clean data from NaN's
        files = files.dropna()
        #Create a list of time-steps
        self.time_list = files['time'].tolist()
        #Create a list of all files' names
        self.list_vtu = files['vtu'].tolist()
        #Format files' names
        self.list_vtu = [i.replace('.pvtu', '.0000.vtu') for i in self.list_vtu]

        if last == None:
            self.list_vtu = self.list_vtu[first::interval]
            self.time_list = self.time_list[first::interval]
        else:
            self.list_vtu = self.list_vtu[first:last:interval]
            self.time_list = self.time_list[first:last:interval]

        #Read VTU data
        N_vtu = len(self.list_vtu)
        pbar = tqdm(N_vtu, desc="Reading VTU files")

        p = Pool(processes=processors)
        numbers = np.arange(len(self.list_vtu)).tolist()
        p.map(self.reader, numbers)
        
        