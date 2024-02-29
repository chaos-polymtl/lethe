from tqdm import tqdm
import shutil

def create_copy(self, prefix):
    # List of paths among read data
    read_files_path_list = [self.pvd_datasets[x].path for x in range(len(self.pvd_datasets))]

    # Write new vtu and pvd files to store modified data.
    # IMPORTANT!!!! If this parameter is empty, that is,  "", data will be written over original vtu and pvd files.
    with open(f'{self.path_output}/{self.pvd_name}') as pvd_in:
        with open(f'{self.path_output}/{prefix}{self.pvd_name}', 'w') as pvd_out:
            for line in pvd_in:
                
                # If line refers to a dataset
                if "vtu" in line:

                    # For all read files
                    for path in read_files_path_list:

                        # If line matches one of the files
                        if path in line:
                            
                            # If vtu is in list_vtu
                            if line.split('file="')[1].split('"/>')[0] in self.list_vtu:
                                line = line.replace('file="', f'file="{prefix}')
                                pvd_out.write(line)
                            read_files_path_list.remove(path)
                            pass
                
                # Write config lines
                else:
                    pvd_out.write(line)

    # Make a copy of VTU files
    n_vtu = len(self.list_vtu)
    pbar = tqdm(total = n_vtu, desc="Writing modified VTU and PVD files")
    new_list_vtu = []
    for i in range(len(self.list_vtu)):
        # Copy file
        shutil.copy2(f'{self.path_output}/{self.list_vtu[i]}', f'{self.path_output}/{prefix}{self.list_vtu[i]}')

        # Append to list of names of VTU files
        new_list_vtu.append(f'{prefix}{self.list_vtu[i]}')
        pbar.update(1)
    
    
    # Fix name of PVD and PVTU files
    self.pvd_name = prefix + self.pvd_name
    self.list_vtu = new_list_vtu