###########################################################
# File : plot_dissipation_rate.py
#---------------------------------------------------------
# Plots the dissipation rate from the kinetic energy
# and from the enstrophy and compares them with a
# solution
###########################################################


import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

font = {#'family' : 'normal',
        'weight' : 'normal',
        'size'   : 13}

plt.rc('font', **font)
colors=['#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02']

prefix_ref="/reference/wang_2013.dat"
fname_ref=os.path.dirname(os.path.realpath(__file__))+prefix_ref

t_ref,ref=np.loadtxt(fname_ref,skiprows=1,unpack=True)


fig = plt.figure(facecolor='white')

ax = fig.add_subplot(111)
plt.plot(t_ref,ref,'-',label="Reference DNS - Wang 2013",color="black",lw=2.)
plt.xlabel('Time [s]')
plt.ylabel('Energy dissipation [W]')
plt.legend(loc=4)

plt.tight_layout()
plt.savefig("reference_comparison.png",dpi=300)
plt.savefig("reference_comparison.pdf")
plt.show()
