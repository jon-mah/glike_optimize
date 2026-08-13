#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks=1

#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=20gb

python simulate_and_fit_ARG.py --OPTIMIZER CMA_ES ../opt_results_CMA_ES/
