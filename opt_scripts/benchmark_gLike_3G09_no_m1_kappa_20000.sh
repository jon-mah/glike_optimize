#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks=1

#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=0
#SBATCH --array=1-10
#SBATCH --output=gL_kappa_20000_%A_%a.log

# SLURM_ARRAY_TASK_ID=1

NUM_TREE=25
NUM_TREE_PREFIX="ntree_${NUM_TREE}_rep_${SLURM_ARRAY_TASK_ID}"

python simulate_and_fit_ARG.py --OPTIMIZER gLike ../opt_results_glike_3G09_no_m1_kappa_20000/$NUM_TREE_PREFIX --NUM_TREES $NUM_TREE --MODEL '3G09_no_m' --SEED 1 --KAPPA 20000
