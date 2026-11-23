#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=bench_gLike
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=30G
#SBATCH --output=bench_true.log

python benchmark_ARG_trees.py --gLike_params ../Analysis
