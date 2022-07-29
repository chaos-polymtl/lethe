import os

with open('case_index.txt') as f:
    mixers = f.readlines()

path = os.getcwd()

for mixer in mixers:
    mixer = mixer.replace('\n','')
    this_path = path + '/' + mixer
    os.chdir(this_path)
    os.system('cp ../launch_mixer.sh .')
    os.system('sbatch -J ' + mixer + ' launch_mixer.sh')
    os.system('cd ../')
