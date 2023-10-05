#############################################################################
"""
Postprocessing code for the capillary wave example.
This code compares the evolution of the relative amplitude of the wave for
different time-step values and the analytical solution from Prosperetti [1].

*** Before running this script, make sure to extract simulation results with
capillary-wave-postprocess.py ***
"""
#############################################################################

#############################################################################
'''Importing Libraries'''
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import itertools
#############################################################################

#############################################################################
# Run script: python3 path_to_capillary_wave-combined.py path_to_analytical_solution_csv_file time-step_multiplier1 time-step_multiplier2...

# Check the number of input arguments
if len(sys.argv)<4:
    print("********************************************************************************\n"
          "Incorrect number of arguments\n"
          "Run script with: \n"
          "\t python3 path_to_capillary_wave-combined.py path_to_analytical_solution_csv_file time-step_multiplier1 time-step_multiplier2...\n"
          "********************************************************************************\n")
    exit(1)

#############################################################################
# Time-step multipliers to plot
time_step_multiplier = sys.argv[2:]

#----------------------------------
# Read quantities of interest
#----------------------------------
# Read Lethe results
i=0
for tms in time_step_multiplier:
    tms=tms.replace('.', '_')
    filename = f"lethe-TSM-{tms}.csv"
    lethe_values = pd.read_csv(filename)
    exec(f"time_list_{i} = np.array(lethe_values['t'])")
    exec(f"relative_amplitude_{i} = np.array(lethe_values['a/a0'])")
    i += 1

# Read verification data from Prosperetti [1]
verification_file_name = sys.argv[1]
analytical_values = pd.read_csv(verification_file_name)
analytical_time = np.array(analytical_values['t'])
analytical_solution = np.array(analytical_values['a/a0'])

#############################################################################
#----------------------------------
# Plot 
#---------------------------------
plt.rcParams['font.size'] = '18'
marker = itertools.cycle(('s', 'x', 'o', '*', 'd'))
fig0 = plt.figure(figsize=(12, 8))
ax0 = fig0.add_subplot(111)
for j in range(0,i):
    tsm_j = time_step_multiplier[j]
    exec(f"plt.plot(time_list_{j}, relative_amplitude_{j}, marker=next(marker), mfc='none', linestyle='', label=f'Lethe — $\Delta t={tsm_j}\Delta t_\sigma$')")
plt.plot(analytical_time, analytical_solution, "--k", linewidth=2, label=f'Prosperetti (1981)')
plt.plot([0,analytical_time[-1]], [1,1], "--", color="tab:gray", linewidth=1)
plt.plot([0,analytical_time[-1]], [-1,-1], "--", color="tab:gray", linewidth=1)
plt.xlim([0,analytical_time[-1]])
plt.xlabel('Time')
plt.ylabel('Relative amplitude')
plt.ticklabel_format(axis='both', style='sci', scilimits=(4,-9))
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('TSM_comparison_figure.png',dpi=300)

#############################################################################
#----------------------------------
# References
#----------------------------------

# [1] A. Prosperetti, “Motion of two superposed viscous fluids,” Phys. Fluids, vol. 24, no. 7, pp. 1217–1223, Jul. 1981, doi: 10.1063/1.863522.