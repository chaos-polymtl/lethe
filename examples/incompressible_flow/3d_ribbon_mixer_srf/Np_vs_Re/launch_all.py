import os

# Read folder names
with open('case_index.txt') as f:
    mixers = f.readlines()

path = os.getcwd()

# Submit a new job for each case
for mixer in mixers:
    mixer = mixer.replace('\n','')
    case_path = path + '/' + mixer
    os.chdir(case_path)
    os.system('sbatch -J ' + mixer + ' launch.sh')      # submit job and name it by its folder's name
    os.system('cd ../')
