#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks=1

#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=20gb
#SBATCH --array=1-20
#SBATCH --output=gLike_4A21_%A_%a.log

# SLURM_ARRAY_TASK_ID=1
NUM_TREE=10
NUM_TREE_PREFIX="ntree_${NUM_TREE}_rep_${SLURM_ARRAY_TASK_ID}"

python simulate_and_fit_ARG.py --OPTIMIZER glike ../opt_results_glike_4A21/$NUM_TREE_PREFIX --NUM_TREES $NUM_TREE --MODEL '4A21' --SEED 1
