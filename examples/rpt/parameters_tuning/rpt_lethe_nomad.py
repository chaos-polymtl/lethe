"""
Name   : rpt_lethe_nomad.py
Author : Audrey Collard-Daigneault
Date   : 17-06-2021
Desc   : This code is a meta Python code to run the rpt_3d lethe application with the NOMAD software for
         blackbox optimization
"""

import os
import sys

# Path to files/application
path = "./"  # Current directory

# Get new calculated parameters
tmpfile = sys.argv[1]  # Temporary file for new values
with open(tmpfile, "r") as values:
    new_values = values.read().split()

# Parameter filenames
initial_lethe_parameter_file = "rpt_parameters.prm"
tmp_lethe_parameter_file = "/tmp/rpt_new_parameters_" + tmpfile.split("/")[2] + ".prm"

# List of parameters to tune
parameters_to_tune = ["dead time", "activity", "attenuation coefficient reactor"]
n_parameters_to_tune = len(parameters_to_tune)

# Search to modify .prm with new calculated tuned parameters
with open(path + initial_lethe_parameter_file, "r") as file:  # Read .prm file with split lines
    filestr = file.read().split("\n")

for line in range(0, len(filestr)):
    if "=" in filestr[line]:
        last_value = filestr[line].split("=")[1]  # get the last value of the line

        for i in range(0, n_parameters_to_tune):
            if parameters_to_tune[i] in filestr[line]:
                filestr[line] = filestr[line].replace(last_value, " " + new_values[i])

with open(tmp_lethe_parameter_file, "w") as file:  # Write .prm with new parameters
    file.write("\n".join(filestr))

# Call rpt_3d executable
os.system("rpt_3d " + tmp_lethe_parameter_file)

# Delete temporary parameter files
os.remove(tmp_lethe_parameter_file)
