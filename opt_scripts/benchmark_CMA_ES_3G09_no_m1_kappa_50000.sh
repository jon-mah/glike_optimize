#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks=1

#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=120gb
#SBATCH --array=1-10
#SBATCH --output=CMA_kappa_50000_%A_%a.log

# SLURM_ARRAY_TASK_ID=1

NUM_TREE=25
NUM_TREE_PREFIX="ntree_${NUM_TREE}_rep_${SLURM_ARRAY_TASK_ID}"

python simulate_and_fit_ARG.py --OPTIMIZER CMA_ES ../opt_results_CMA_ES_3G09_no_m1_kappa_50000/$NUM_TREE_PREFIX --NUM_TREES $NUM_TREE --MODEL '3G09_no_m' --SEED 1 --KAPPA 50000
