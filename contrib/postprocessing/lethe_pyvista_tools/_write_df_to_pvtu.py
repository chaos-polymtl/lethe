from tqdm import tqdm


# Write modifications on each df to PVTU files
def write_df_to_pvtu(self, prefix = "mod_"):
    """
    Writes .pvd and .pvtu files from data stored in self.df.
    The files are written in self.output_path.

    Parameter:

    :param prefix = "mod_"           -> String with prefix of the written files.
    By default, "mod_" is added in front of the regular files.
    """

    # List of paths among read data
    read_files_path_list = [self.pvd_datasets[x].path for x in range(len(self.pvd_datasets))]

    # Write modified PVD to match new PVTU files
    if len(prefix) > 0:
        with open(f'{self.path_output}/{self.pvd_name}') as pvd_in:
            with open(f'{self.path_output}/{prefix}{self.pvd_name}', 'w') as pvd_out:
                for line in pvd_in:
                    
                    # If line refers to a dataset
                    if "pvtu" in line:

                        # For all read files
                        for path in read_files_path_list:

                            # If line matches one of the files
                            if path in line:
                                                               
                                # If pvtu is in list_pvtu
                                if line.split('file="')[1].split('"/>')[0] in self.list_pvtu:
                                    line = line.replace('file="', f'file="{prefix}')
                                    pvd_out.write(line)
                                read_files_path_list.remove(path)
                                pass
                    
                    # Write config lines
                    else:
                        pvd_out.write(line)
    
    
    if self.df_available:
        # Write modified PVTU file
        N_pvtu = len(self.df)
        pbar = tqdm(total = N_pvtu, desc="Writing new PVTU and PVD files")
        for i in range(len(self.df)):
            self.df[i].save(f'{self.path_output}/{prefix}{self.list_vtu[i]}', binary = False)
            pbar.update(1)

        print(f"Modified .pvtu and .pvd files with prefix {prefix} successfully written")
    
    else:
        print(f"No df available for writing. Try to use read_lethe_to_pyvista first")