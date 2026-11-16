#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks=1

#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=0
#SBATCH --array=1-10
#SBATCH --output=CMA_10000_true_%A_%a.log

# SLURM_ARRAY_TASK_ID=1

NUM_TREE=25
NUM_TREE_PREFIX="ntree_${NUM_TREE}_rep_${SLURM_ARRAY_TASK_ID}"

python simulate_and_fit_ARG.py --OPTIMIZER CMA_ES ../opt_results_CMA_ES_3G09_no_m1_kappa_10000_true/$NUM_TREE_PREFIX --NUM_TREES $NUM_TREE --MODEL '3G09_no_m' --SEED 1 --KAPPA 10000 --TRUE_START --START_SEED $SLURM_ARRAY_TASK_ID
