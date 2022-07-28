#!/bin/bash
#SBATCH --time=5:00:00
#SBATCH --account=rrg-blaisbru
#SBATCH --ntasks-per-node=40
#SBATCH --nodes=1
#SBATCH --mem=180G
#SBATCH --output=%x-%j.out

source $HOME/.dealii
srun  gls_navier_stokes_3d ribbon_gls.prm
