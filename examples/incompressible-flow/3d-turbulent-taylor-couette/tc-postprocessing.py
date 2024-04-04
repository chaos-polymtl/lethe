# -*- coding: utf-8 -*-

import numpy as np
import argparse
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import os 
from tc_functions import *


#Plot font and colors
font = {'weight' : 'normal',
        'size'   : 13}

plt.rc('font', **font)
colors=['#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02']

#Parameter of the simulation 
class parameter(): 
    re = 1.0          #Outer radius
    ri = 0.5          #Inner radius
    kappa = 0.5       #Radius ratio (ri/re)
    omega = 1.0       #Angular velocity (rad/s)
    d = 0.5           #Width of the annulus (re - ri)
    alpha = 2*np.pi   #Aspect ratio (L/d)
    nu = 6.25e-5      #Kinematic viscosity
    rho = 1.0         #Fluid density 
    Re = 4000         #Reynolds number

prm = parameter() 

### Simulations data ###

#Load enstrophy file
parser = argparse.ArgumentParser(description='Arguments for calculation of the dissipation rate')
parser.add_argument("-ens", "--enstrophy", type=str, help="Name of the input file for the enstrophy", required=True)
parser.add_argument("-l", "--legend", type=str, help="Legend label", required=False)
args, leftovers=parser.parse_known_args()
fname_ens=args.enstrophy
leg = args.legend


t_1, e = enstrophy(prm,fname_ens)
t_p3, e_ref_p3, t_p4, e_ref_p4, t_p5, e_ref_p5 = enstrophy_ref()
 

### Plots section ###
plt.close('all')


# Enstrophy #
plt.figure(2)

if (leg==None):
   plt.plot(t_1,e,'k', label = "Lethe",lw=2)
else:
   plt.plot(t_1,e,'k', label = leg,lw=2)
plt.plot(t_p3,e_ref_p3 ,label="Wang P3",lw=2, linestyle='--', color=colors[0])
plt.plot(t_p4,e_ref_p4 ,label="Wang P4",lw=2, linestyle='--', color=colors[1])
plt.plot(t_p5,e_ref_p5 ,label="Wang P5",lw=2, linestyle='--', color=colors[2])
plt.xlabel("Time [-]")
plt.ylabel("Enstrophy [-]")
plt.legend(loc = 'lower right')
# Zoom #
ax_zoom = plt.axes([0.18, 0.65, 0.2, 0.2])  
ax_zoom.yaxis.set_major_locator((MaxNLocator(integer=True)))
ax_zoom.plot(t_1,e,'k',lw=2)
plt.plot(t_p3,e_ref_p3 , linestyle='--',lw=2, color=colors[0])
plt.plot(t_p4,e_ref_p4 , linestyle='--',lw=2, color=colors[1])
plt.plot(t_p5,e_ref_p5 , linestyle='--',lw=2, color=colors[2])
ax_zoom.set_xlim(16, 40)
ax_zoom.set_ylim(12, 18)
ax_zoom.xaxis.set_major_locator((MaxNLocator(nbins=4,integer=True)))
ax_zoom.yaxis.set_major_locator((MaxNLocator(nbins=4,integer=True)))
if (leg==None):
  plt.savefig('enstrophy_comparison.png', dpi=300)
else:
  plt.savefig(leg.replace(" ","")+".png", dpi=300)
plt.show()