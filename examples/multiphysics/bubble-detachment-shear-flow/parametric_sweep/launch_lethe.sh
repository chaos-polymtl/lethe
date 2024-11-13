#!/bin/bash
#SBATCH --account=acount-name
#SBATCH --ntasks-per-node=64
#SBATCH --nodes=1 
#SBATCH --mem=150G
#SBATCH --time=24:00:00
#SBATCH --job-name=sweep
#SBATCH --mail-type=BEGIN 
#SBATCH --mail-type=END 
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=your-email@your-mailbox

source $HOME/.dealii
srun $HOME/lethe/inst/bin/lethe-fluid bubble-detachment-shear-flow.prm
