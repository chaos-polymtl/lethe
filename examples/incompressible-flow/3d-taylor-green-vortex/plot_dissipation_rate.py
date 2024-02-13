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

parser = argparse.ArgumentParser(description='Arguments for calculation of the dissipation rate')
parser.add_argument("-ke", "--kinetic_rate", type=str, help="Name of the input file for kinetic energy dissipation rate", required=True)
parser.add_argument("-ens", "--enstrophy", type=str, help="Name of the input file for the enstrophy", required=True)

parser.add_argument("-v", "--viscosity", type=float, help="viscosity", required=True)
parser.add_argument("-z", "--zoom", type=bool, help="viscosity", required=False)
args, leftovers=parser.parse_known_args()

viscosity=args.viscosity
fname_ke=args.kinetic_rate
fname_ens=args.enstrophy
zoom=args.zoom

prefix_ref="/reference/wang_2013.dat"
fname_ref=os.path.dirname(os.path.realpath(__file__))+prefix_ref

t_ike,ike=np.loadtxt(fname_ke,skiprows=1,unpack=True)
t_ens,ens=np.loadtxt(fname_ens,skiprows=1,unpack=True)
t_ref,ref=np.loadtxt(fname_ref,skiprows=1,unpack=True)

ens=ens*viscosity*2 

fig = plt.figure(facecolor='white')

ax = fig.add_subplot(111)
plt.plot(t_ike,ike,label="Kinetic energy dissipation",color=colors[1],lw=2.)
plt.plot(t_ens,ens,label="Enstrophy energy dissipation",color=colors[0],lw=2.)
plt.plot(t_ref,ref,'--',label="Reference",color="black",lw=2.)
plt.xlabel('Time [s]')
plt.ylabel('Energy dissipation [W]')
plt.legend(loc=4)

if (zoom):
    #this is an inset axes over the main axes
    inset_axes = inset_axes(ax, 
                    width="25%", # width = 30% of parent_bbox
                    height=1.4, # height : 1 inch
                    loc=1)
    # the main axes is subplot(111) by default
    plt.plot(t_ike,ike,label="Kinetic energy dissipation",color=colors[1],lw=2)
    plt.plot(t_ens,ens,label="Enstrophy energy dissipation",color=colors[0],lw=2)
    plt.plot(t_ref,ref,'--',label="Reference",color="black",lw=2)
    plt.axis([8.0,10.0,0.011,0.013])

plt.tight_layout()
    
plt.savefig("dissipation_comparison.png",dpi=300)
plt.savefig("dissipation_comparison.pdf")
plt.show()
