# SPDX-FileCopyrightText: Copyright (c) 2024 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#############################################################################
"""
Postprocessing code for thermocapillary migration of a bubble example

"""
#############################################################################

'''Importing Libraries'''
import numpy as np
import sys
import matplotlib.pyplot as plt
import argparse
#############################################################################

radius = 0.25
gamma_prime = -1

lambda_ratio = 1
beta_ratio = 1

rho = 1

nu_ext = 1
nu_int = nu_ext*lambda_ratio

mu_ext = nu_ext*rho
mu_int = nu_int*rho


#Take case path as argument and store it
parser = argparse.ArgumentParser(description='Arguments for calculation of the dissipation rate')
parser.add_argument("-i", "--input",nargs='+', help="Name of the input files", required=True)
parser.add_argument("-l", "--labels",nargs='+', help="Labels in the legend, must be equal to the length of the input files", required=True)
args, leftovers=parser.parse_known_args()

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
for i in range(len(args.input)):
    output_dir = sys.argv[i+2]
    filename = output_dir + "/barycenter_information.dat"
    
    t,x,y, z,vx,vy, vz=np.loadtxt(filename,skiprows=1,unpack=True)
    ax1.plot(t, vx , '-', lw=2, label="Lethe, " + args.labels[i])


t_analytical = np.linspace(0.0,3.0, 25)

vx_analytical = -2*radius*gamma_prime/(mu_ext*(2+3*lambda_ratio)*(2+beta_ratio))*np.ones(len(t_analytical))

ax1.plot(t_analytical, vx_analytical, 'ok',label="Analytical velocity")
# plt.yscale('log')
# plt.xscale('log')

plt.xlim([0, 3])
ax1.set_ylabel(r'Migration velocity')
ax1.set_xlabel(r'$t$')
ax1.legend(loc="upper left")
ax1.legend(loc=4)
fig1.savefig(f'./bubble-migration-velocity.png',dpi=300)
plt.show()


