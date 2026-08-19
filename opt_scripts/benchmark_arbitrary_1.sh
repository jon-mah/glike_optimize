#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks=1

#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=20gb
#SBATCH --array=1-5

NUM_TREE=1
NUM_TREE_PREFIX="ntree_${NUM_TREE}_rep_${SLURM_ARRAY_TASK_ID}"

python simulate_and_fit_ARG.py --OPTIMIZER glike ../opt_results_glike_arbitrary_rep/$NUM_TREE_PREFIX --NUM_TREES $NUM_TREE
