# SPDX-FileCopyrightText: Copyright (c) 2024 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#############################################################################
"""
Postprocessing code for the Rayleigh-Plateau example.
This code calculates an average breakup length for the different cases and
compares them with simulation results from Denner et al. [1]
"""
#############################################################################

#############################################################################
'''Importing Libraries'''
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
sys.path.append("$LETHE_PATH/contrib/postprocessing/")
from lethe_pyvista_tools import *
#############################################################################

#############################################################################
# Run script: python3 path_to_rayleigh-plateau-compare.py
#             path_to_reference_solution_csv_file delta_value1 delta_value2...

# Check the number of input arguments
if len(sys.argv)<4:
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
          "Incorrect number of arguments\n"
          "Run script with: \n"
          "\t python3 path_to_rayleigh-plateau-compare.py path_to_reference_solution_csv_file delta_value1 delta_value2...\n "
          "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    exit(1)

#############################################################################
# Excitation amplitudes to plot
delta_values = sys.argv[2:]

# Initialize lists
Lb_list = []
delta_list = []
no_breakup_delta_list = []

# Jet radius
r_jet = 1.145e-3

#----------------------------------
# Read quantities of interest
#----------------------------------
# Read Lethe results
for delta in delta_values:
    delta = float(delta)
    delta_string = (f"{delta:.2f}").replace('.', '_')
    filename = f"lethe-delta{delta_string}.csv"
    lethe_values = pd.read_csv(filename)
    breakup_lengths = np.array(lethe_values['Lb'])
    if len(breakup_lengths) > 2:
        delta_list.append(delta)
        breakup_lengths = np.delete(breakup_lengths, [0,1]) # We remove the first 2 breakups
        Lb_list.append(np.mean(breakup_lengths))
    else:
        no_breakup_delta_list.append(delta)


# Get dimensionless breakup lengths
length_list = [Lb / r_jet for Lb in Lb_list]

# Read reference data from Denner et al. [1]
reference_data_file = sys.argv[1]
reference_data = pd.read_csv(reference_data_file)
reference_delta_list = np.array(reference_data['delta0'])
reference_length_list = np.array(reference_data['Lb_r0'])

#############################################################################
#----------------------------------
# Plot
#----------------------------------
plt.rcParams['font.size'] = '18'
fig0 = plt.figure(figsize=(12, 8))
if len(delta_list) > 0:
    plt.plot(delta_list, length_list, "sr", label="Lethe (2024)")
for i, delta_value in enumerate(no_breakup_delta_list):
    if i == 0:
        plt.axvline(x=delta_value, linestyle='--', linewidth=2, color="r",
                    label="no breakup - Lethe (2024)")
    else:
        plt.axvline(x=delta_value, linestyle='--', linewidth=2, color="r")
plt.plot(reference_delta_list, reference_length_list, "xb", markeredgewidth=2,
         label="Denner et al. (2022)")
plt.xlabel(r"Excitation amplitude ($\delta_0$)")
plt.ylabel(r"$L_b/R_{jet}$")
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('rayleigh-plateau_comparison_figure.png',dpi=300)
plt.show()

#############################################################################
#----------------------------------
# References
#----------------------------------

# [1] F. Denner, F. Evrard, A. A. Castrejón-Pita, J. R. Castrejón-Pita, and B. van Wachem, “Reversal and Inversion of Capillary Jet Breakup at Large Excitation Amplitudes,” Flow Turbul. Combust., vol. 108, no. 3, pp. 843–863, Mar. 2022, doi: 10.1007/s10494-021-00291-w.
