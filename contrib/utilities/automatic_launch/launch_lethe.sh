#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --account=$yourgroupaccount
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=32G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=$your.email.adress@email.provider
#SBATCH --output=%x-%j.out

source $HOME/.dealii
srun $HOME/lethe/inst/bin/gls_navier_stokes_2d cylinder.prm