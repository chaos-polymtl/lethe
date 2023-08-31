from tqdm import tqdm

# Sort all data given reference array 
def sort_by_array(self, reference_array_name = "ID"):
    '''
    Sorts all self.df according to a reference array:

    Parameter:

    reference_array_name = "ID" -> String with name of reference array.
    "ID" is used as default for particles, but any other 1D array can be 
    used.
    '''
    
    pbar = tqdm(total = len(self.time_list), desc = f"Sorting dataframe by {reference_array_name}")
    
    if self.df_available:
        for i in range(len(self.time_list)):
            self.df[i].points = self.df[i].points[self.df[i][reference_array_name].argsort()]
            for name in self.df[0].array_names:
                self.df[i][name] = self.df[i][name][self.df[i][reference_array_name].argsort()]
            pbar.update(1)
    
    else:
        for i in range(len(self.list_vtu)):
            df = self.get_df(i)
            df.points = df.points[df[reference_array_name].argsort()]
            for name in df.array_names:
                df[name] = df[name][df[reference_array_name].argsort()]
            
            df.save(f'{self.path_output}/{self.list_vtu[i]}')
            pbar.update(1)