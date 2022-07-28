# -*- coding: utf-8 -*-
"""
Script used to prepare multiple simulations subdirectories and modify .prm to replace a 
-------------------
Preparation: put an initial .prm file along with needed .msh and .sh in a directory.
-------------------
 * key : key identifying the string to replace using jinja with the value of the parameters
 * values : list of values to use for each of the parameters

Output :
 * subdirectories with simulation files needed and modified .prm
 * txt file with the name of each subdirectory
"""

import jinja2
import shutil
import os
import numpy as np

#system characteristics
D=0.27                                      # impeller diameter
N= 10 / 2 / np.pi                           # angular velocity of the agitator
Re = np.logspace(-1,np.log10(100),25)       # Re values (different cases)
nu = N * D * D / Re                         # kinematic vicosity values


folder_prefix="mixer_"
values=nu

template_folder="template"
template_prm_file="ribbon_gls.prm"

templateLoader = jinja2.FileSystemLoader(searchpath="./")
templateEnv = jinja2.Environment(loader=templateLoader)
template = templateEnv.get_template(template_prm_file)


case_index = open('case_index.txt', "w")


for val in values:

    copy_folder = folder_prefix + str(val)
    if os.path.exists(copy_folder) and os.path.isdir(copy_folder):
        shutil.rmtree(copy_folder)
    shutil.copytree(template_folder, copy_folder)

    # Copy folder    
    outputText = template.render(N=val)
    with open(copy_folder +"/"+ template_prm_file, 'w') as f:
        f.write(outputText)
    
    case_index.write(copy_folder+"\n")

