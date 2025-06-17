#!/bin/bash -x
#SBATCH --cluster=hpda2
#SBATCH --partition=hpda2_compute_gpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:4
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --mail-user=bgeihe@uni-koeln.de
#SBATCH --mail-type=all
#SBATCH --job-name=single_node
#SBATCH --output=stdout.%j
#SBATCH --error=stderr.%j
#SBATCH --time=00:30:00

source profile

mpiexec -n $SLURM_NTASKS $JL --threads=1 --project=. run.jl

