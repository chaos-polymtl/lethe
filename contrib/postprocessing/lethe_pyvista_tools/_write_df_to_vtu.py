from tqdm import tqdm


# Write modifications on each df to VTU files
def write_df_to_vtu(self, prefix = "mod_"):
    """
    Writes .pvd and .vtu files from data stored in self.df.
    The files are written in self.output_path.

    Parameter:

    :param prefix = "mod_"           -> String with prefix of the written files.
    By default, "mod_" is added in front of the regular files.
    """

    # List of paths among read data
    read_files_path_list = [self.pvd_datasets[x].path for x in range(len(self.pvd_datasets))]

    # Write modified PVD to match new VTU files
    with open(f'{self.path_output}/{self.pvd_name}') as pvd_in:
        with open(f'{self.path_output}/{prefix}{self.pvd_name}', 'w') as pvd_out:
            for line in pvd_in:
                
                # If line refers to a dataset
                if "vtu" in line:

                    # For all read files
                    for path in read_files_path_list:

                        # If line matches one of the files
                        if path in line:
                            line = line.replace('.pvtu', '.0000.vtu')
                            
                            # If vtu is in list_vtu
                            if line.split('file="')[1].split('"/>')[0] in self.list_vtu:
                                line = line.replace('file="', f'file="{prefix}')
                                pvd_out.write(line)
                            read_files_path_list.remove(path)
                            pass
                
                # Write config lines
                else:
                    pvd_out.write(line)
    
    
    if self.df_available:
        # Write modified VTU file
        N_vtu = len(self.df)
        pbar = tqdm(total = N_vtu, desc="Writing new VTU and PVD files")
        for i in range(len(self.df)):
            self.df[i].save(f'{self.path_output}/{prefix}{self.list_vtu[i]}')
            pbar.update(1)

        print(f"Modified .vtu and .pvd files with prefix {prefix} successfully written")
    
    else:
        print(f"No df available for writing. Try to use read_lethe_to_pyvista first")