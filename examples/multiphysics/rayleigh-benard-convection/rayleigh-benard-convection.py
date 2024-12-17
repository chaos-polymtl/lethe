# SPDX-FileCopyrightText: Copyright (c) 2022-2024 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#############################################################################
"""
Postprocessing code for rayleigh-benard-convection example

"""
#############################################################################

'''Importing Libraries'''
import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv
import argparse
import os
import re

#############################################################################

#Functions to read the data from the reference files
def get_last_file(directory):

    #Check if the directory exists
    if not os.path.isdir(directory):
        print(f"Directory '{directory}' not found")
        return None

    # Create the regular expression pattern for the data files
    pattern = re.compile(r"rayleigh-benard_convection\.(\d+)\.00000\.vtu")

    # Get a list of the files
    files = [
        f for f in os.listdir(directory) 
        if os.path.isfile(os.path.join(directory, f)) and pattern.match(f)
    ]
    
    # Check if any files were found
    if not files:
        print(f"No .vtu files found in '{directory}' matching the regex pattern.")
        return None  

    # Sort files by modification time
    sorted_files = sorted(files, key=lambda x: int(x.split('.')[1]))
    
    # Return the most recently modified file with the given extension
    return sorted_files[-1]

#############################################################################


#Data from Ouertatani, Cheikh, Beya, Lili (2008) for Ra= 10^4
x_ref_Ouertatani = [0.000800640512409928,0.00960768614891914,0.0168134507606085,0.0256204963971177,0.0344275420336269,0.044035228182546,0.0544435548438751,0.0640512409927942,0.0760608486789432,0.0880704563650921,0.100080064051241,0.11369095276221,0.127301841473179,0.140912730184147,0.155324259407526,0.170536429143315,0.185748598879103,0.203362690152122,0.21857485988791,0.236188951160929,0.254603682946357,0.273819055244195,0.292233787029624,0.311449159327462,0.3306645316253,0.349879903923139,0.371497197758207,0.392313851080865,0.412329863891113,0.433146517213771,0.453963170536429,0.476381104883907,0.497998398718975,0.519615692554043,0.540432345876701,0.562850280224179,0.582866293034428,0.602882305844676,0.624499599679744,0.645316253002402,0.66453162530024,0.685348278622898,0.703763010408327,0.722978382706165,0.740592473979183,0.759007205764612,0.77662129703763,0.792634107285829,0.809447558046437,0.825460368294636,0.841473178542834,0.856685348278623,0.870296236989592,0.88390712570056,0.897518014411529,0.911128903122498,0.919935948759007,0.933546837469976,0.944755804643715,0.953562850280224,0.963971176941553,0.973578863090472,0.981585268214572,0.990392313851081,0.99679743795036]
y_ref_Ouertatani = [-0.00088809946714032, 0.00532859680284192, 0.0133214920071048, 0.022202486678508, 0.0328596802841918, 0.0426287744227354, 0.0506216696269982, 0.0612788632326821, 0.0746003552397869, 0.0870337477797513, 0.0976909413854352, 0.11101243339254, 0.124333925399645, 0.13854351687389, 0.153641207815275, 0.168738898756661, 0.184724689165187, 0.200710479573712, 0.217584369449378, 0.235346358792185, 0.252220248667851, 0.271758436944938, 0.291296625222025, 0.309946714031972, 0.329484902309059, 0.348134991119005, 0.369449378330373, 0.390763765541741, 0.412078152753108, 0.432504440497336, 0.454706927175844, 0.475133214920071, 0.497335701598579, 0.518650088809947, 0.539076376554174, 0.560390763765542, 0.58259325044405, 0.603907637655417, 0.623445825932504, 0.644760213143872, 0.6651865008881, 0.683836589698046, 0.704262877442274, 0.723801065719361, 0.740674955595027, 0.759325044404973, 0.77708703374778, 0.792184724689165, 0.809946714031972, 0.825932504440497, 0.841030195381883, 0.857904085257549, 0.871225577264654, 0.884547069271759, 0.897868561278863, 0.911190053285968, 0.921847246891652, 0.935168738898757, 0.94404973357016, 0.953818827708703, 0.965364120781528, 0.974245115452931, 0.984014209591474, 0.990230905861457, 1]
x_vel_Ouertatani = [-0.000647249190938504, -0.0213592233009709, -0.0433656957928803, -0.0666666666666667, -0.0899676375404531, -0.111974110032362, -0.133333333333333, -0.154045307443366, -0.171521035598706, -0.187702265372168, -0.202588996763754, -0.214886731391586, -0.226537216828479, -0.236893203883495, -0.244012944983819, -0.248543689320388, -0.251132686084142, -0.250485436893204, -0.249838187702265, -0.245307443365696, -0.238834951456311, -0.229126213592233, -0.218770226537217, -0.205177993527508, -0.188996763754045, -0.170226537216829, -0.150809061488673, -0.128802588996764, -0.106148867313916, -0.0815533980582525, -0.0563106796116505, -0.0304207119741101, -0.00258899676375407, 0.0233009708737864, 0.0498381877022653, 0.0757281553398058, 0.0996763754045307, 0.123624595469256, 0.146278317152104, 0.165695792880259, 0.184466019417476, 0.200647249190938, 0.214886731391586, 0.226537216828479, 0.236893203883495, 0.244660194174757, 0.248543689320388, 0.252427184466019, 0.251779935275081, 0.250485436893204, 0.245954692556634, 0.238834951456311, 0.229773462783171, 0.218770226537217, 0.205825242718447, 0.191585760517799, 0.176051779935275, 0.157281553398058, 0.137864077669903, 0.117799352750809, 0.0957928802588997, 0.0744336569579288, 0.0504854368932039, 0.0271844660194175, 0.00194174757281551]
y_vel_Ouertatani = [-0.00039564787339263496,0.023343224530168183,0.04787339268051438,0.07240356083086058,0.09851632047477749, 0.12383778437190907, 0.1452027695351138, 0.1673590504451039, 0.18793273986152326, 0.20613254203758657, 0.22116716122650848, 0.2354104846686449, 0.24727992087042538, 0.25598417408506435, 0.25994065281899115, 0.26310583580613256, 0.2623145400593472, 0.25994065281899115, 0.2536102868447082, 0.24490603363006924, 0.23303659742828886, 0.22116716122650848, 0.2053412462908012, 0.18872403560830864, 0.1713155291790307, 0.152324431256182, 0.13254203758654803, 0.11196834817012857,0.09060336300692384,0.06844708209693373,0.047082096933729,0.02571711177052427,0.0011869436201780714, -0.02096933728981204, -0.042334322453016826, -0.06369930761622156, -0.08664688427299705, -0.1072205736894164, -0.12858555885262118, -0.14757665677546983, -0.1681503461918892, -0.1847675568743818, -0.2029673590504451, -0.2172106824925816, -0.23066271018793275, -0.24253214638971315, -0.2536102868447082, -0.25914935707220577, -0.2623145400593472, -0.26468842729970327, -0.2623145400593472, -0.257566765578635, -0.2504451038575668, -0.237784371909001, -0.22670623145400595, -0.21088031651829872, -0.19347181008902078, -0.17289812067260138, -0.152324431256182, -0.1301681503461919, -0.10642927794263107,-0.08189910979228487,-0.056577645895153295,-0.03046488625123639,-0.0027695351137487223]

#Take case path as argument and store it
parser = argparse.ArgumentParser(description='Arguments for the validation of the rayleigh-benard-convection example')
parser.add_argument("--validate", action="store_true", help="Launches the script in validation mode. This will log the content of the graph and prevent the display of figures", default=False)
parser.add_argument("-f", "--folder", type=str, help="Path to the output folder. This is the folder that contains the results of the simulation (.vtu, .pvtu, .dat and .pvd files)", required=True)
args, leftovers=parser.parse_known_args()
output_dir =args.folder

#Load the data from the result folder
filename = output_dir + "/"+ get_last_file(output_dir)
#sim = pv.read(filename)
sim = pv.read("/home/oreste/work/lethe/lethe/examples/multiphysics/rayleigh-benard-convection/output/rayleigh-benard_convection.193600.pvtu")
sim.set_active_vectors("velocity")

#Extract the data on the middle line in x and y direction
a = [-0.5,0,0]
b = [0.5,0,0]
sampled_data_x = sim.sample_over_line(a, b, resolution=100)

a = [0,-0.5,0]
b = [0,0.5,0]
sampled_data_y = sim.sample_over_line(a, b, resolution=100)

#Get the position and the velocity of the simulation
x = sampled_data_x["Distance"]
y = sampled_data_y["Distance"]
u = sampled_data_y["velocity"][:,0]
v = sampled_data_x["velocity"][:,1]


sim = pv.read("/home/oreste/work/lethe/lethe/examples/multiphysics/rayleigh-benard-convection/output/rayleigh-benard_convection.20000.pvtu")
sim.set_active_vectors("velocity")

#Extract the data on the middle line in x and y direction
a = [-0.5,0,0]
b = [0.5,0,0]
sampled_data_x = sim.sample_over_line(a, b, resolution=100)

a = [0,-0.5,0]
b = [0,0.5,0]
sampled_data_y = sim.sample_over_line(a, b, resolution=100)

#Get the position and the velocity of the simulation
x2 = sampled_data_x["Distance"]
y2 = sampled_data_y["Distance"]
u2 = sampled_data_y["velocity"][:,0]
v2 = sampled_data_x["velocity"][:,1]

print(u-u2)
print(v-v2)

exit()


#Normalize the velocity by the characteristic velocity u_char = sqrt(g*L*beta*(T_h-T_c))
u /= np.sqrt(10*1*0.71*10)
v /= np.sqrt(10*1*0.71*10)

fig0 = plt.figure()
ax0 = fig0.add_subplot(111)
ax0.plot(-u, y, lw=2, label="Lethe")
ax0.plot(x_vel_Ouertatani,y_ref_Ouertatani, 's',ms=8,mfc="None",markeredgecolor="black",label="Ouertatani 2008")

ax0.set_ylabel(r'$y$')
ax0.set_xlabel(r'$u$')
ax0.set_xlim(-0.5,0.5)
ax0.set_ylim(0,1)
ax0.set_aspect('equal')
ax0.legend(loc="upper left")
if (args.validate):
  solution = np.column_stack((u, y))
  np.savetxt("solution-rayleigh_10k-uy.dat",solution, header="u y")
  fig0.savefig(f'./solution-rayleigh_10k-uy.pdf')
else:
  fig0.savefig(f'./solution-rayleigh_10k-uy.png',dpi=300)
  plt.show()

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.plot(x, -v, lw=2, label="Lethe")
ax1.plot(x_ref_Ouertatani, y_vel_Ouertatani, 's',ms=8,mfc="None",markeredgecolor="black",label="Ouertatani 2008")

ax1.set_ylabel(r'$v$')
ax1.set_xlabel(r'$x$')
ax1.set_xlim(0,1)
ax1.set_ylim(-0.5,0.5)
ax1.set_aspect('equal')
ax1.legend(loc="upper right")
if (args.validate):
  solution = np.column_stack((x, v))
  np.savetxt("solution-rayleigh_10k-xv.dat",solution, header="x v")
  fig1.savefig(f'./solution-rayleigh_10k-xv.pdf')

else:
  fig1.savefig(f'./solution-rayleigh_10k-xv.png',dpi=300)
  plt.show()

