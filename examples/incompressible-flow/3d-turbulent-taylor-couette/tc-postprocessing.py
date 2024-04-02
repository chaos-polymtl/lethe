# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import os 
from tc_functions import *



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

#Enstrophy

t_1, e = enstrophy(prm,'enstrophy.dat')
t_p3, e_ref_p3, t_p4, e_ref_p4, t_p5, e_ref_p5 = enstrophy_ref()
 

### Plots section ###
plt.close('all')


# Enstrophy #
plt.figure(2)
plt.plot(t_1,e,'k', label = "Lethe")
plt.plot(t_p3,e_ref_p3,'g',label="Wang P3", linestyle='--')
plt.plot(t_p4,e_ref_p4,color='orange',label="Wang P4", linestyle=':')
plt.plot(t_p5,e_ref_p5 ,'r',label="Wang P5", linestyle='-.')
plt.xlabel("Time [-]")
plt.ylabel("Enstrophy [-]")
plt.legend(loc = 'lower right')
# Zoom #
ax_zoom = plt.axes([0.18, 0.65, 0.2, 0.2])  
ax_zoom.yaxis.set_major_locator((MaxNLocator(integer=True)))
ax_zoom.plot(t_1,e,'k')
ax_zoom.plot(t_p3,e_ref_p3, 'g', linestyle='--')
ax_zoom.plot(t_p4,e_ref_p4, color='orange', linestyle=':')
ax_zoom.plot(t_p5,e_ref_p5,'r',linestyle='-.')
ax_zoom.set_xlim(16, 40)
ax_zoom.set_ylim(12, 18)
ax_zoom.xaxis.set_major_locator((MaxNLocator(nbins=4,integer=True)))
ax_zoom.yaxis.set_major_locator((MaxNLocator(nbins=4,integer=True)))
plt.savefig('enstrophy_comparison.png', dpi=300)