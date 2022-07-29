#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --account=$yourgroupaccount
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40 
#SBATCH --output=%x-%j.out

source $HOME/.dealii
srun  $HOME/lethe/inst/bin/gls_navier_stokes_3d ribbon_gls.prm
