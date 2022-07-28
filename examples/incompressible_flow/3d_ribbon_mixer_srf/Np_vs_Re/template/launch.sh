#!/bin/bash
#SBATCH --time=23:30:00
#SBATCH --account=rrg-blaisbru
#SBATCH --nodes=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=blais.bruno@gmail.com
#SBATCH --ntasks-per-node=40 
#SBATCH --output=%x-%j.out

source $HOME/.dealii
srun  gls_navier_stokes_3d ribbon_gls.prm
