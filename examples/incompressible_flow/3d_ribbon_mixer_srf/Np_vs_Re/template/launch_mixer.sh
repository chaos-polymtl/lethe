#!/bin/bash
#SBATCH --account=
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --mem=90G
#SBATCH --output=%x-%j.out

source $HOME/.dealii
srun  gls_navier_stokes_3d ribbon_gls.prm
