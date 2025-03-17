"""
Summary: Script to launch different cases (as folders) of the same problem 
with different levels of refinement and polynomial degrees for the velocity and the pressure.
"""
import os

PATH = os.getcwd()
PRM_FILE = 'MMS_analysis_2D_steady_auto.prm'
LETHE_EXEC = 'mpirun -np 6 lethe-fluid'

for root, directories, files in os.walk(PATH):
    if PRM_FILE in files and root != PATH:
                os.chdir(root)
                os.system(f'{LETHE_EXEC} {PRM_FILE}' + '> output.txt')