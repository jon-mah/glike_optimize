#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks=1

#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=50G
#SBATCH --array=1-10
#SBATCH --output=CMA_fixed_bug_%A_%a.log

# SLURM_ARRAY_TASK_ID=1

NUM_TREE=20
NUM_TREE_PREFIX="ntree_${NUM_TREE}_rep_${SLURM_ARRAY_TASK_ID}"

python simulate_and_fit_ARG.py --OPTIMIZER CMA_ES ../opt_results_CMA_ES_3G09_no_m1_kappa_10000_bug/$NUM_TREE_PREFIX --NUM_TREES $NUM_TREE --MODEL '3G09_no_m' --SEED 1 --KAPPA 10000 --trees-dir /project2/jitang_1167/projects/glike_sglike_data/glike_on_3g09_model/trees_reps20-trees20_msprime1.0.0_fix-all-mig-rates --rep $SLURM_ARRAY_TASK_ID
# python simulate_and_fit_ARG.py --OPTIMIZER CMA_ES ../opt_results_CMA_ES_3G09_no_m1_kappa_10000/$NUM_TREE_PREFIX --NUM_TREES $NUM_TREE --MODEL '3G09_no_m' --SEED 1 --KAPPA 10000 --rep $SLURM_ARRAY_TASK_ID
