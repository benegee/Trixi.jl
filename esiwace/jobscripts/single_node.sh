#!/bin/bash -x
#SBATCH --cluster=hpda2
#SBATCH --partition=hpda2_compute_gpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:4
#SBATCH --mail-user=<mail_addr>
#SBATCH --mail-type=all
#SBATCH --export=NONE
#SBATCH --output=stdout.%j
#SBATCH --error=stderr.%j
#SBATCH --time=00:30:00

source profile

srun $JL --threads=1 --project=. run.jl
