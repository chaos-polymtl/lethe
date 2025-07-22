# SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

######################################################################
# Import Libraries

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import argparse
from matplotlib.colors import LogNorm

# Set plot parameters

plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['figure.figsize'] = (10,8)
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['lines.markeredgewidth'] = 2
plt.rcParams['legend.columnspacing'] = 2
plt.rcParams['legend.handlelength'] = 1
plt.rcParams['legend.handletextpad'] = 0.2
plt.rcParams['legend.frameon'] = True
plt.rcParams['legend.fancybox'] = False
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['xtick.major.size'] = 5
plt.rcParams['ytick.major.size'] = 5
plt.rcParams['ytick.major.width'] = 2
plt.rcParams['font.size'] = '20'
plt.rcParams['font.family']='DejaVu Serif'
plt.rcParams['font.serif']='cm'
plt.rcParams['savefig.bbox']='tight'
plt.rcParams.update({
    'text.usetex': True,
    'text.latex.preamble': r'\usepackage{amsfonts}'
})

import sys
sys.path.append("$LETHE_PATH/contrib/postprocessing/")
from lethe_pyvista_tools import *

######################################################################


parser = argparse.ArgumentParser(description='Arguments for the post-processing of the Pseudo-2d gas solid fluidized bed example')
parser.add_argument("-f", "--folder", type=str, help="Folder path. This folder is the folder which contains the .prm file.", required=True)
args, leftovers=parser.parse_known_args()

# Simulation folder
folder=args.folder

# Get data
data = pd.read_csv(folder + "velocity_disturbance.dat")

# pi*alpha^2*dp/h where alpha = (3/(4pi))^(1/3) ; pi*alpha^2 ~ 6/5
d = np.logspace(-4,1,100)
upper_bound = 6./5. * d

# Plot data
fig, ax = plt.subplots()
norm = LogNorm(vmin=0.01, vmax=500)
sc = ax.scatter(data['dp/h'], data['|u_h-u_inf|/u_inf'], c=data['Rep'], cmap='jet', norm=norm, alpha=0.8, edgecolors='k')
ax.set_xscale('log')
ax.set_yscale('log')
cb = fig.colorbar(sc, ax=ax)
cb.set_label(r'$Re_p$')

plt.loglog(d,upper_bound,'k--',label=r"$\pi\left(3/(4\pi)\right)^{2/3}d_p/h$")
plt.xlim(1e-4, 1e1)
# plt.ylim(1e-4,1e1)
plt.xlabel(r"$d_p/h$")
plt.ylabel(r"$\frac{|u_h-u_{\infty}|}{u_\infty}$")
plt.grid()
plt.legend()
plt.savefig('velocity_disturbance')
plt.show() 
