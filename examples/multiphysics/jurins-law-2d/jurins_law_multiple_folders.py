# SPDX-FileCopyrightText: Copyright (c) 2023-2024 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

"""
Postprocessing code for Jurin's law example
This code extracts the difference in height between the meniscus
and the fluid on the side.

Quick user guide: 

This python script is used to plot the difference in height between the meniscus
and the fluid on the side for different cases of capillary rise.

How to use this script?

1- Create a first directory (we will call it outputs but feel free to name it as you like).
2- For each simulation, create a new directory in outputs : output1, output2, output3,...
3- Run the simulations and store the results in their correct directories
4- Execute the script with the path of outputs as argument
5- Enjoy your plots. Each curve will be given the name of the directory containing the data required to plot it. In this example, I will get 3 curves. with names output1, output2 and output3.
"""
# -------------------------------------------
# Modules
# -------------------------------------------

from postprocessing_jurins_law_dimensioned import get_deltaH, analytical_solution
import os
import sys
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from matplotlib.lines import Line2D
import scienceplots
import numpy as np
from natsort import os_sorted

plt.style.use(['science','ieee'])

#For controlling font sized globally
SMALL_SIZE = 8
MEDIUM_SIZE = 8
BIGGER_SIZE = 15

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', labelsize=12)             # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


rootdir =sys.argv[1]

folder_name_list = []
root, dirs, files = next(os.walk(rootdir, topdown=True))

for dir in dirs:
    folder_name_list.append(str(root + "/" + dir))
    
folder_name_list = os_sorted(folder_name_list)
dirs = os_sorted(dirs)

class parametres():
    g=9.810
    phase_limit=0.0
    mu = 3.00e-2
    sigma = 7.3e-2
    rho_l = 2.00e3
    r = 1e-3
    
    
prm=parametres()

pparams=dict(xlabel=r'$t \text{[s]}$', ylabel=r'$\Delta H \text{[mm]}$')

with plt.style.context(['science','ieee']):
      

   fig = plt.figure()
   ax = fig.add_subplot(111)
   #For continuous color spectrum
   color = iter(cm.seismic(np.linspace(0, 1, len(folder_name_list))))

   for i in range(len(folder_name_list)):
       c = next(color)
       deltaH, time_list = get_deltaH(folder_name_list[i],prm)
       label_loop = dirs[i]
       if (label_loop=="90"):
          ax.plot(time_list, deltaH,lw=1,label=label_loop +"°", linestyle='solid',color='k')
       else:
          angle = float(label_loop[:2])
          a_solution = analytical_solution(prm,angle)
          plt.axhline(y=a_solution,color='k',ls='dotted')
          ax.plot(time_list, deltaH,lw=1,label=label_loop +"°", linestyle='solid',color=c)
       
   handles, labels = plt.gca().get_legend_handles_labels()
   line = Line2D([0], [0], label='Analytical', color='k',ls='dotted')
   handles.append(line)
   box = ax.get_position()
   ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
   ax.set_ylim([-4,4])
   fig.legend(loc='outside center right',frameon = True,edgecolor='k',prop={'size': MEDIUM_SIZE},ncol=1, fancybox=False, bbox_to_anchor=(1.13, 0.5),handles=handles)
   ax.set(**pparams)
   fig.savefig('height_differences.pdf',format="pdf",dpi=500)
   fig.savefig('height_differences.png',format="png",dpi=500)
