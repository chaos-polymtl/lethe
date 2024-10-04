# SPDX-FileCopyrightText: Copyright (c) 2024 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

import sympy as sym
import argparse
import numpy as np

parser = argparse.ArgumentParser(description='Arguments for calculation of the dissipation rate')
parser.add_argument("-f", "--folder", type=str, help="Folder path", required=True)
parser.add_argument("-i", "--input", type=str, help="Name of the input file", required=True)
parser.add_argument("-o", "--output", type=str, default="ke_rate.dat",  help="Name of the output file", required=False)
args, leftovers=parser.parse_known_args()

folder=args.folder
filename=args.input
outname=args.output

data=np.loadtxt(folder + filename,skiprows=1)

# Sybmolic position x
x = sym.Symbol('x')

# Data in x: time step
positions=data[:,0]

# Data in y: kinetic energy
y=data[:,1]

size = len(positions)

# Function that calculates derivative at node i using a Lagrange polynomial 
# of degree deg and a stencil from left_index to right_index. It works for 
# equidistant (fixed time step) and non-equidistant nodes (adaptive time
# step).
def calculate_derivative(deg, i, left_index, right_index):
    v = [1]*(deg + 1) 
    v_derivative = [0]*(deg + 1)
    result = 0
    derivative = 0
    count = 0
    for k in range(left_index, right_index + 1):
        for j in range(left_index, right_index + 1):
            if k !=j:
                # l_k(x)
                v[count] *= (x-positions[j])/(positions[k]-positions[j])

        #Calculate derivative of l_k(x)
        v_derivative[count] = sym.diff(v[count],x)
        
        #Calculate L(x) as a linear combination of y_i*l_k(x)    
        result += y[k]*v[count]
        derivative += y[k]*v_derivative[count]
        count += 1

    # Evaluate the Lagrange polynomial at the desired x position
    return derivative.evalf(subs={x:positions[i]}) 

final_ke_rate=[]
for i in range(0,size):
    kerate=0.0
    if i==0:
        kerate = calculate_derivative(1, i, i, i+1)
    elif i==size-1:
        kerate = calculate_derivative(1, i, i-1, i)
    elif (i==1) or (i==(size-2)):
        kerate = calculate_derivative(2, i, i-1, i+1)
    else:
        kerate = calculate_derivative(4, i, i-2, i+2)

    final_ke_rate.append(kerate* -1)

# Save data with the time step and the kinetic energy rate in two columns
DataOut = np.column_stack((positions,final_ke_rate))
np.savetxt(folder + outname,DataOut)
