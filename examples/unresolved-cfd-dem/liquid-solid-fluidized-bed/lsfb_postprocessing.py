
# To use the lethe_pyvista_tools, you need to have python 3
# installed in your computer

# The modules necessary to run lethe pyvista tools are:
# Pandas: pip install pandas
# PyVista: pip install pyvista

# To use lethe_pyvista_tools, add /contrib/postprocessing to your PythonPATH 
# adding 
# ``export PYTHONPATH="${PYTHONPATH}:$LETHE_PATH/contrib/postprocessing"``
# to your ~/.bashrc, or append the path to the /contrib/postprocessing
# folder in Lethe to sys

import sys
# Path to the module
path_to_module = '$LETHE_PATH/contrib/postprocessing/'
sys.path.append(path_to_module)

# or even put the "lethe_pyvista_tools.py" file inside
# the same directory as your python script and procceed as follows

# This line imports all lethe_pyvista_tools functionalities
from lethe_pyvista_tools import *


#############################################################################
'''Parameters'''

figure = False #True if plot
plot_P_t = False
write_to_excel = True
#############################################################################

#############################################################################
'''Importing Libraries'''
from math import pi
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn import linear_model
from tqdm import tqdm

import os
#############################################################################

#############################################################################
'''Simulation properties'''

#Take case path as argument and store at case_path
case_path = sys.argv[1]
saveFigDir = case_path


fluid = lethe_pyvista_tools(case_path = case_path, prm_file_name = "liquid-solid-fluidized-bed.prm", pvd_name = "cfd_dem.pvd", step = 100)

particle = lethe_pyvista_tools(case_path = case_path, prm_file_name = "liquid-solid-fluidized-bed.prm", pvd_name = "cfd_dem_particles.pvd", step = 100)

dp = particle.prm_dict['diameter']
rp = dp/2
Vp = 4/3*pi*rp**3
rhop = particle.prm_dict['density particles']
Np = particle.prm_dict['number'][1]
inlet_velocity = fluid.prm_dict['u']


g = abs(particle.prm_dict['gx'])       #m/s^2
Vnp = Np * 4/3 * pi * rp**3

rhol = fluid.prm_dict['density']  #kg/m^3
nu = fluid.prm_dict['kinematic viscosity'] #m^2/s
mu = rhol*nu    #Pa
Db = 0.1        #m
Rb = Db/2       #m
Area = pi*Rb**2 #m^2
Hb = 1.1

#############################################################################

#############################################################################
'''Functions'''
def Analytical(Np, rp, rhol, rhop, g, Area):
    Mp = Np*4/3*pi*rp**3*rhop
    delta_p = Mp*(rhop - rhol)*g/(rhop * Area)
    return delta_p
#############################################################################

#Define eps_list and voidfraction_list to append value for each time-step
eps_list = []
voidfraction_list = []
total_deltaP = []

#Create a list with all x values for the pressure takes
x_pressure_takes = np.arange(-0.45, Hb/2 + 0.01, 0.01)

deltaP_analytical = Analytical(Np, rp, rhol, rhop, g, Area)

if os.path.isdir(case_path + '/P_x') == False:
    os.mkdir(case_path + '/P_x')

pbar = tqdm(total = len(fluid.time_list)-1, desc="Processing data")
for i in range(len(fluid.time_list)):
    #Define "df" as last time step df
    df = fluid.get_df(i) 

    #Take the first slice of the domain at pos0
    pos0 = [-0.45, 0, 0]
    slice0 = df.slice(normal=[1, 0, 0], origin = pos0)
    p0 = slice0.integrate_data()['pressure']/slice0.area

    #Create empty lists to fill with x values and pressure as function of x values
    p_x = []
    x_list = []

    #Loop through z values above z0
    for j in range(len(x_pressure_takes)):
        x_list.append(x_pressure_takes[j]-pos0[0])
        slice_x = df.slice(normal=[1, 0, 0], origin = [x_pressure_takes[j], 0, 0])
        p_x.append((p0 - slice_x.integrate_data()['pressure']/slice_x.area)[0]*rhol)

    #Store the total pressure drop
    total_deltaP.append(p_x[-1])

    #Apply least-squares to find the porosity
    regr = linear_model.LinearRegression() #fit_intercept = 0
    x_list = np.array(x_list).reshape(-1, 1)
    p_x = np.array(p_x).reshape(-1, 1)
    model = regr.fit(x_list, p_x)
    r_sq = model.score(x_list, p_x)

    #Find linear portion of the graph and adjust a linear function to it
    x = x_list
    y = p_x

    for j in range(len(y)-1):
        if y[j] >= deltaP_analytical*0.9:
                x = x[:j]
                y = y[:j]
                break
    
    model = regr.fit(x, y)

    #Calculate bed voidage by the slope of the pressure drop curve
    eps = 1-(model.coef_[0][0]/((rhop-rhol)*g))
    eps_list.append(eps)


    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)

    fig1.suptitle(f'Time = {fluid.time_list[i]} s')
    ax1.plot(x_list, p_x, 'ok', ms = 5, label ='Simulation')
    ax1.plot(x_list, np.repeat(deltaP_analytical, len(x_list)), '.r', ms = 2, label = 'Analytical')
    ax1.plot(x, y, '.g')
    ax1.plot(x, model.predict(x), '-b')
    ax1.set_ylabel(r'$-\Delta p \/\ [Pa]$')
    ax1.set_xlabel(r'$Height \/\ [m]$')
    ax1.set_xlim(0, 1.02)
    ax1.set_ylim(0, deltaP_analytical*1.30)
    ax1.legend()
    ax1.annotate(r'$\varepsilon (-dp/dz) = {:1.2}$'.format(eps), (x[round(len(x)/2)], y[round(len(y)/2)]), xytext=(0.65, 0.4), textcoords='axes fraction', arrowprops=dict(facecolor='black', shrink=0.04), fontsize=14, horizontalalignment='right', verticalalignment='top')
    fig1.savefig(f'{saveFigDir}/P_x/P_x-{i}.png')
    plt.close(fig1)

    pbar.update(1)

#Export the total pressure results as a function of time
csv = pd.DataFrame([fluid.time_list, total_deltaP], index=['time', 'deltaP']).transpose()
csv.to_csv(f'{saveFigDir}/deltaP_t.csv')

print(f'Average - delta p = {np.mean(total_deltaP)} Pa')


#Plot void fraction with time
fig0, ax0 = plt.subplots(1, 1)

ax0.plot(fluid.time_list, eps_list, 'ok', ms = 5)
ax0.legend()
ax0.grid()
#ax0.set_xlim(0, 35)
ax0.set_ylim(0, 1)
ax0.set_ylabel(r'$Bed \/\ voidage \/\ [-]$')
ax0.set_xlabel(r'$Time \/\ [s]$')
fig0.savefig(f'{saveFigDir}/eps_t.png')
plt.close(fig0)

eps_ave = np.average(eps_list[-10:])
eps_std = np.std(eps_list[-10:])

print('Average void fraction among the last 10 time steps: ' + str(eps_ave))


