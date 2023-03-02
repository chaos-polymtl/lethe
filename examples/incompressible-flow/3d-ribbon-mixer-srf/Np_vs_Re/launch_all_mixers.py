import os

PATH = os.getcwd()

# User input
CASE_PREFIX = 'mixer_nu_'
PRM_FILE = 'ribbon-gls.prm'
SHELL_FILE = 'launch-mixer.sh'

for root, directories, files in os.walk(PATH):

    if CASE_PREFIX in root:

        os.chdir(root)

        case_name = root.split('/')[-1]
        os.system(f'sbatch -J {case_name} {SHELL_FILE}')