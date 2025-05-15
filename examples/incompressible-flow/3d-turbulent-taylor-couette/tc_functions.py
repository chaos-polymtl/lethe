# SPDX-FileCopyrightText: Copyright (c) 2024 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

import numpy as np 
import os


def enstrophy (prm,file_name):
    """Function that retrieves and returns, in the form of arrays, the time and dimensionless kinetic energy values obtained during simulation.  

    Inputs:
        - prm : Object class parameter () 
            - ri = Radius of the inner cylinder
            - rho = Fluid density
            - omega = Angular velocity of the inner cylinder
            - d = Annulus width (re-ri)
           
    Outputs (same order as below):
        - t1 = Vector (np.array) of simulation times
        - E1 = Vector (np.array) of all dimensionless kinetic energy values
    """ 
    
    #Parameter variables
    rho = prm.rho
    omega = prm.omega
    ri = prm.ri
    d = prm.d
    U = omega*ri
    
    #Dimensionless factor
    coeff = (2*d*d)/(rho*U*U)
    
    #Recovering data
    file = file_name 
    path = os.path.join(os.path.dirname(__file__), file)

    #Loading data
    t1,E1 = np.loadtxt(path, skiprows =1 , unpack = True)

    #Adimensionnal enstrophy 
    E1 = E1 * coeff
    
    return np.array(t1[10:]),np.array(E1[10:])
    
def enstrophy_ref(): 
    """Function that retrieves and returns, in the form of arrays, the time and kinetic energy values of references.  

        
            
    Outputs (same order as below):
        - t4 = Vector (np.array) of reference times (P4)
        - E4 = Vector (np.array) of all reference dimensionless enstrophy values (P4)
        - t5 = Vector (np.array) of reference times (P5)
        - E5 = Vector (np.array) of all reference dimensionless enstrophy values (P5)
    """ 
    
    #Recovering data
    file3 = 'enstrophy_wang_p3.tsv'
    file4 = 'enstrophy_wang_p4.tsv'
    file5 = 'enstrophy_wang_p5.tsv'
    dir = 'references'
    path3 = os.path.join(os.path.dirname(__file__),dir,file3)
    path4 = os.path.join(os.path.dirname(__file__),dir,file4)
    path5 = os.path.join(os.path.dirname(__file__),dir,file5)
    
    #Loading data
    t3,E3 = np.loadtxt(path3, skiprows = 1, unpack = True) 
    t4,E4 = np.loadtxt(path4, skiprows = 1, unpack = True)
    t5,E5 = np.loadtxt(path5, skiprows = 1, unpack = True)

    return np.array(t3),np.array(E3), np.array(t4),np.array(E4),np.array(t5),np.array(E5)



