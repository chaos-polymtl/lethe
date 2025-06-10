# SPDX-FileCopyrightText: Copyright (c) 2022-2025 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#############################################################################
"""
Postprocessing code for rising-bubble example

"""
#############################################################################

'''Importing Libraries'''
import numpy as np
import sys
import matplotlib.pyplot as plt
import pyvista as pv
import argparse
import os
import pandas as pd
#############################################################################

'''Plot formating'''

from cycler import cycler

colors=['#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02']

plt.rcParams['axes.prop_cycle'] = cycler(color = colors)
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['figure.figsize'] = (10,8)
plt.rcParams['lines.linewidth'] = 4
plt.rcParams['lines.markersize'] = '11'
plt.rcParams['markers.fillstyle'] = "none"
plt.rcParams['lines.markeredgewidth'] = 2
plt.rcParams['legend.columnspacing'] = 2
plt.rcParams['legend.handlelength'] = 2.8
plt.rcParams['legend.handletextpad'] = 0.2
plt.rcParams['legend.frameon'] = True
plt.rcParams['legend.fancybox'] = False
plt.rcParams['legend.fontsize'] = '18'
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['xtick.major.size'] = 5
plt.rcParams['ytick.major.size'] = 5
plt.rcParams['ytick.major.width'] = 2
# plt.rcParams['axes.margins'] = 15
plt.rcParams['font.size'] = '25'
plt.rcParams['font.family']='DejaVu Serif'
plt.rcParams['font.serif']='cm'
plt.rcParams['savefig.bbox']='tight'

plt.rcParams.update({
    'text.usetex': True,
    'text.latex.preamble': r'\usepackage{amsfonts}'
})

#############################################################################

parser = argparse.ArgumentParser(description='Arguments for the validation of the 2d lid-driven cavity')
parser.add_argument("--validate", action="store_true", help="Launches the script in validation mode. This will log the content of the graph and prevent the display of figures", default=False)
parser.add_argument("-s", "--sharp", type=str, help="Path to the output folder for the sharpening results. This is the folder that contains the results of the simulation (.vtu, .pvtu, .dat and .pvd files)", required=True)
parser.add_argument("-g", "--geo", type=str, help="Path to the output folder for the geometric results. This is the folder that contains the results of the simulation (.vtu, .pvtu, .dat and .pvd files)", required=True)
parser.add_argument("-a", "--alge", type=str, help="Path to the output folder for the algebraic results. This is the folder that contains the results of the simulation (.vtu, .pvtu, .dat and .pvd files)", required=True)

parser.add_argument("-c", "--case", type=int, choices=[1, 2], help="Test case number. Can be either 1 or 2", required=True)
args, leftovers=parser.parse_known_args()
case_number = args.case

output_dir_sharp =args.sharp
filename_barycenter_sharp = output_dir_sharp + "/vof_barycenter_information.dat"
filename_mass_sharp = output_dir_sharp + "/mass_conservation_information.dat"

output_dir_geo =args.geo
filename_barycenter_geo = output_dir_geo + "/vof_barycenter_information.dat"
filename_mass_geo = output_dir_geo + "/mass_conservation_information.dat"

output_dir_alge =args.alge
filename_barycenter_alge = output_dir_alge + "/vof_barycenter_information.dat"
filename_mass_alge = output_dir_alge + "/mass_conservation_information.dat"

# Get the data from the references for test case 1
#Data from Zahedi, Kronbichler and Kreiss (2012)
x_ref_ZKR = [0.047, 0.139, 0.232 ,0.324 ,0.416 ,0.508 ,0.601 ,0.693 ,0.785 ,0.878 ,0.97 ,1.062 ,1.154 ,1.247 ,1.339 ,1.431 ,1.524 ,1.616 ,1.708 ,1.8 ,1.893 ,1.985 ,2.077 ,2.17 ,2.262 ,2.354 ,2.446 ,2.539 ,2.631 ,2.723 ,2.816 ,2.908 ,2.979]
y_ref_ZKR = [0.501 ,0.505 ,0.513 ,0.525 ,0.54 ,0.557 ,0.576 ,0.597 ,0.618 ,0.641 ,0.663 ,0.685 ,0.707 ,0.729 ,0.75 ,0.77 ,0.791 ,0.811 ,0.83 ,0.849 ,0.868 ,0.886 ,0.904 ,0.922 ,0.94 ,0.958 ,0.976 ,0.993 ,1.011 ,1.029 ,1.047 ,1.064 ,1.078]
x_vel_ZKR = [0.001 ,0.03 ,0.056 ,0.081 ,0.106 ,0.136 ,0.17 ,0.203 ,0.237 ,0.271 ,0.304 ,0.338 ,0.372 ,0.41 ,0.452 ,0.494 ,0.544 ,0.607 ,0.687 ,0.779 ,0.872 ,0.964 ,1.056 ,1.149 ,1.241 ,1.333 ,1.426 ,1.518 ,1.61 ,1.702 ,1.795 ,1.887 ,1.979 ,2.072 ,2.164 ,2.256 ,2.348 ,2.441 ,2.533 ,2.625 ,2.718 ,2.81 ,2.902 ,2.978]
y_vel_ZKR = [0.004 ,0.017 ,0.029 ,0.041 ,0.053 ,0.066 ,0.081 ,0.095 ,0.109 ,0.122 ,0.135 ,0.147 ,0.158 ,0.17 ,0.182 ,0.194 ,0.205 ,0.218 ,0.229 ,0.237 ,0.24 ,0.241 ,0.239 ,0.236 ,0.231 ,0.227 ,0.222 ,0.217 ,0.213 ,0.208 ,0.205 ,0.201 ,0.198 ,0.196 ,0.194 ,0.192 ,0.192 ,0.191 ,0.191 ,0.192 ,0.192 ,0.193 ,0.193 ,0.194]
df_contour_ZKR = pd.read_csv("reference_contour/case1_contour_ZKR.csv", header=0)

#Data Hysing, S., Turek, S., Kuzmin, D., Parolini, N., Burman, E., Ganesan, S., & Tobiska, L. (2009). Quantitative benchmark computations of two‐dimensional bubble dynamics. International Journal for Numerical Methods in Fluids, 60(11), 1259-1288.
x_ref_H = [0.24476 , 0.49894 , 0.75132 , 0.99538 , 1.2466  , 1.49797 , 1.75014 , 1.99467 , 2.24561 , 2.49821 ,2.74912 ]
y_ref_H = [0.514646, 0.554469, 0.608872, 0.670196, 0.728744, 0.785978, 0.838403, 0.888369, 0.936482, 0.984353, 1.032372]
x_vel_H = [0.24535, 0.49537, 0.74679   ,0.99191   ,1.2424,   1.49542   ,1.74728   ,1.99275   ,2.24447   ,2.49706,   2.74796]
y_vel_H = [0.114461, 0.196094, 0.23631, 0.241153, 0.231292, 0.218334, 0.206795, 0.197974, 0.193075, 0.191356, 0.192096]


# Get the data from the references for test case 2
# Data TP2D in Hysing, S., Turek, S., Kuzmin, D., Parolini, N., Burman, E., Ganesan, S., & Tobiska, L. (2009). Quantitative benchmark computations of two‐dimensional bubble dynamics. International Journal for Numerical Methods in Fluids, 60(11), 1259-1288.
x_ref_TP2D = [0.002, 0.246, 0.493, 0.736, 0.988, 1.233, 1.485, 1.729, 1.975, 2.226, 2.471, 2.717, 3.000]
y_ref_TP2D = [0.501, 0.517, 0.563, 0.623, 0.685, 0.741, 0.796, 0.853, 0.909, 0.968, 1.024, 1.077, 1.138]
x_vel_TP2D = [0.000, 0.245, 0.498, 0.748, 0.992, 1.243, 1.496, 1.748, 1.994, 2.245, 2.495, 2.746, 2.999]
y_vel_TP2D = [0.001, 0.139, 0.230, 0.253, 0.242, 0.229, 0.226, 0.232, 0.243, 0.232, 0.220, 0.217, 0.217]
df_contour_TP2D = pd.read_csv("reference_contour/case2_contour_TP2D.csv", header=0)

# Data FreeLIFE in Hysing, S., Turek, S., Kuzmin, D., Parolini, N., Burman, E., Ganesan, S., & Tobiska, L. (2009). Quantitative benchmark computations of two‐dimensional bubble dynamics. International Journal for Numerical Methods in Fluids, 60(11), 1259-1288.
x_ref_FreeLIFE = [0.005, 0.246, 0.498, 0.750, 1.000, 1.246, 1.499, 1.750, 2.000, 2.247, 2.496, 2.748, 3.000]
y_ref_FreeLIFE = [0.501, 0.517, 0.563, 0.625, 0.685, 0.742, 0.797, 0.853, 0.912, 0.969, 1.023, 1.073, 1.125]
x_vel_FreeLIFE = [0.000, 0.244, 0.496, 0.748, 0.998, 1.245, 1.496, 1.747, 1.999, 2.243, 2.496, 2.747, 2.999]
y_vel_FreeLIFE = [0.001, 0.138, 0.229, 0.252, 0.239, 0.227, 0.224, 0.232, 0.244, 0.225, 0.212, 0.208, 0.205]
df_contour_FreeLIFE = pd.read_csv("reference_contour/case2_contour_FreeLIFE.csv", header=0)

# Data MooNMD in Hysing, S., Turek, S., Kuzmin, D., Parolini, N., Burman, E., Ganesan, S., & Tobiska, L. (2009). Quantitative benchmark computations of two‐dimensional bubble dynamics. International Journal for Numerical Methods in Fluids, 60(11), 1259-1288.
x_ref_MooNMD = [0.002, 0.246, 0.493, 0.736, 0.988, 1.233, 1.485, 1.729, 1.975, 2.226, 2.471, 2.717, 3.000]
y_ref_MooNMD = [0.501, 0.517, 0.563, 0.623, 0.685, 0.741, 0.796, 0.853, 0.909, 0.968, 1.024, 1.077, 1.138]
x_vel_MooNMD = [0.246, 0.487, 0.735, 0.986, 1.231, 1.482, 1.728, 1.973, 2.226, 2.471, 2.715, 0.002, 2.999]
y_vel_MooNMD = [0.137, 0.226, 0.250, 0.239, 0.226, 0.222, 0.228, 0.238, 0.231, 0.218, 0.215, 0.001, 0.213]
df_contour_MooNMD = pd.read_csv("reference_contour/case2_contour_MooNMD.csv", header=0)


# Load position and velocity of the barycenter
t_sharp,x_sharp,y_sharp,vx_sharp,vy_sharp=np.loadtxt(filename_barycenter_sharp,skiprows=1,unpack=True)
t_geo,x_geo,y_geo,vx_geo,vy_geo=np.loadtxt(filename_barycenter_geo,skiprows=1,unpack=True)
t_alge,x_alge,y_alge,vx_alge,vy_alge=np.loadtxt(filename_barycenter_alge,skiprows=1,unpack=True)


# Plot the position of the barycenter
fig0 = plt.figure()
ax0 = fig0.add_subplot(111)
ax0.plot(t_sharp, y_sharp, '-', lw=2, label="Projection")
ax0.plot(t_geo, y_geo, '-', lw=2, label="Geometric")
ax0.plot(t_alge, y_alge, '-', lw=2, label="PDE-based")


if case_number == 1 :

  ax0.plot(x_ref_ZKR, y_ref_ZKR, 'ok',label="Zahedi et al. (2012)")
  ax0.plot(x_ref_H, y_ref_H, 'sk',alpha=0.8,label="Hysing et al. (2009)")

elif case_number ==2 :

  ax0.plot(x_ref_MooNMD, y_ref_MooNMD, '^k',alpha=0.8,label="MooNMD, Hysing et al. (2009)")
  ax0.plot(x_ref_FreeLIFE, y_ref_FreeLIFE, 'sk',alpha=0.8,label="FreeLIFE, Hysing et al. (2009)")
  ax0.plot(x_ref_TP2D, y_ref_TP2D, 'ok',alpha=0.8, label="TP2D, Hysing et al. (2009)")

ax0.set_ylabel(r'Barycenter height [L]')
ax0.set_xlabel(r'Rising time [T]')
ax0.legend(loc="upper left")
fig0.savefig(f'./ymean-t-case' + str(case_number) + '.png',dpi=300)
plt.show()

# Plot the velocity
fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.plot(t_sharp, vy_sharp, '-', lw=2, label="Projection")
ax1.plot(t_geo, vy_geo, '-', lw=2, label="Geometric")
ax1.plot(t_alge, vy_alge, '-', lw=2, label="PDE-based")

if case_number == 1 :

  ax1.plot(x_vel_ZKR, y_vel_ZKR, 'ok',label="Zahedi et al. (2012)")
  ax1.plot(x_vel_H, y_vel_H, 'sk',alpha=0.8,label="Hysing et al. (2009)")

elif case_number == 2 :

  ax1.plot(x_vel_MooNMD, y_vel_MooNMD, '^k',alpha=0.8,label="MooNMD, Hysing et al. (2009)")
  ax1.plot(x_vel_FreeLIFE, y_vel_FreeLIFE, 'sk',alpha=0.8,label="FreeLIFE, Hysing et al. (2009)")
  ax1.plot(x_vel_TP2D, y_vel_TP2D, 'ok', alpha=0.8,label="TP2D, Hysing et al. (2009)")
  
ax1.set_ylabel(r'Rise velocity [LT$^{-1}$]')
ax1.set_xlabel(r'Rising time [T]')
ax1.legend(loc=4)
fig1.savefig(f'./bubble-rise-velocity-case' + str(case_number) + '.png',dpi=300)
plt.show()

# Plot the contour of the bubble in the last frame
list_vtu = os.listdir(output_dir_sharp)
list_vtu = [(output_dir_sharp+"/"+x) for x in list_vtu if  ("vtu" in x)]
latest_file = max(list_vtu, key=os.path.getctime)
if(not args.validate): print("Opening file: ", latest_file)
sim = pv.read(latest_file)
sim.set_active_scalars("filtered_phase")

contour_val = np.array([0.5])
contours = sim.contour(contour_val)
x_sharp, y_sharp = contours.points[:, 0], contours.points[:, 1]

list_vtu = os.listdir(output_dir_geo)
list_vtu = [(output_dir_geo+"/"+x) for x in list_vtu if  ("vtu" in x)]
latest_file = max(list_vtu, key=os.path.getctime)
if(not args.validate): print("Opening file: ", latest_file)
sim = pv.read(latest_file)
sim.set_active_scalars("filtered_phase")

contour_val = np.array([0.5])
contours = sim.contour(contour_val)
x_geo, y_geo = contours.points[:, 0], contours.points[:, 1]

list_vtu = os.listdir(output_dir_alge)
list_vtu = [(output_dir_alge+"/"+x) for x in list_vtu if  ("vtu" in x)]
latest_file = max(list_vtu, key=os.path.getctime)
if(not args.validate): print("Opening file: ", latest_file)
sim = pv.read(latest_file)
sim.set_active_scalars("filtered_phase")

contour_val = np.array([0.5])
contours = sim.contour(contour_val)
x_alge, y_alge = contours.points[:, 0], contours.points[:, 1]


x = [x_sharp, x_geo, x_alge]
y = [y_sharp, y_geo, y_alge]

label = ['Projection', 'Geometric', 'PDE-based']
fig_name = ['sharp', 'geo', 'alge']


for i in range(3):
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    
    ax2.scatter(x[i], y[i], s=2, marker="s", color=colors[i], label=label[i], linewidths=1)
    
    if case_number == 1 :

      ax2.plot(df_contour_ZKR['x'], df_contour_ZKR['y'], '-k',
        alpha=0.8,label="Zahedi et al. (2012)", linewidth=1)
      ax2.set_ylim([0.8,1.4])

    elif case_number == 2 :
      ax2.plot(df_contour_MooNMD['x'], df_contour_MooNMD['y'], '-k',
              alpha=0.8,label="MooNMD, Hysing et al. (2009)", linewidth=1)
              
      ax2.plot(df_contour_FreeLIFE['x'], df_contour_FreeLIFE['y'], '--k',
              alpha=0.8,label="FreeLIFE, Hysing et al. (2009)", linewidth=1)
      ax2.plot(df_contour_FreeLIFE['x_1'], df_contour_FreeLIFE['y_1'], '--k',
              alpha=0.8,linewidth=1)
      ax2.plot(df_contour_FreeLIFE['x_2'], df_contour_FreeLIFE['y_2'], '--k',
              alpha=0.8,linewidth=1)
              
      ax2.plot(df_contour_TP2D['x'], df_contour_TP2D['y'], '-.k',
              alpha=0.8,label="TP2D, Hysing et al. (2009)", linewidth=1)
      ax2.plot(df_contour_TP2D['x_1'], df_contour_TP2D['y_1'], '-.k',
              alpha=0.8,linewidth=1)
      ax2.plot(df_contour_TP2D['x_2'], df_contour_TP2D['y_2'], '-.k',
              alpha=0.8,linewidth=1)
      ax2.set_xlim([0,1])
      ax2.set_ylim([0.5,1.4])

    ax2.set_xlabel(r'x [L]')
    ax2.set_ylabel(r'y [L]')
    ax2.legend(markerscale=1, scatterpoints=20)
    ax2.grid( which='major', color='grey', linestyle='--')

    if (args.validate):
      solution = np.column_stack((x, y))
      np.savetxt("solution-contour-case" + str(case_number) + ".dat", solution, header="x y")
      fig2.savefig("bubble-contour-case" + str(case_number) + ".png")
    else:
      fig2.savefig(fig_name[i] + "-bubble-contour-case" + str(case_number) + ".png",dpi=300)
      plt.show()

# Plot the mass conservation for the global and local mass
df_mass_sharp = pd.read_csv(filename_mass_sharp, header=0, sep=r'\s+')
df_mass_geo = pd.read_csv(filename_mass_geo, header=0, sep=r'\s+')
df_mass_alge = pd.read_csv(filename_mass_alge, header=0, sep=r'\s+')

fig3 = plt.figure()
ax3 = fig3.add_subplot(111)
ax3.plot(df_mass_sharp['time'], df_mass_sharp['surface_fluid_1']/df_mass_sharp['surface_fluid_1'].iloc[0], "-", label=r'Projection', linewidth=2)
ax3.plot(df_mass_geo['time'], df_mass_geo['surface_fluid_1']/df_mass_geo['surface_fluid_1'].iloc[0], "-", label=r'Geometric', linewidth=2)
ax3.plot(df_mass_alge['time'], df_mass_alge['surface_fluid_1']/df_mass_alge['surface_fluid_1'].iloc[0], "-", label=r'PDE-based', linewidth=2)

ax3.set_ylabel(r'$V/{V_\mathrm{0}}[-]$')
ax3.set_xlabel(r'Rising time [T]')
ax3.ticklabel_format(useOffset=False, style='plain', axis='y')
ax3.legend(loc="upper left")
fig3.savefig(f'./global-mass-conservation-case' + str(case_number) + '.png',dpi=300)
plt.show()

fig4 = plt.figure()
ax4 = fig4.add_subplot(111)
ax4.plot(df_mass_sharp['time'], df_mass_sharp['geometric_surface_fluid_1']/df_mass_sharp['geometric_surface_fluid_1'].iloc[0], "-", label=r'Projection',linewidth=2)
ax4.plot(df_mass_geo['time'], df_mass_geo['geometric_surface_fluid_1']/df_mass_geo['geometric_surface_fluid_1'].iloc[0], "-", label=r'Geometric',linewidth=2)
ax4.plot(df_mass_alge['time'], df_mass_alge['geometric_surface_fluid_1']/df_mass_alge['geometric_surface_fluid_1'].iloc[0], "-", label=r'PDE-based',linewidth=2)
ax4.set_ylabel(r'$V/{V_\mathrm{0}}[-]$')
ax4.set_xlabel(r'Rising time [T]')
ax4.ticklabel_format(useOffset=False, style='plain', axis='y')
ax4.legend(loc="upper left")
fig4.savefig(f'./geo-mass-conservation-case' + str(case_number) + '.png',dpi=300)
plt.show()



