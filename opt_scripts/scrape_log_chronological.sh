#!/bin/bash
#SBATCH --array=1-10
#SBATCH --time=00:10:00
#SBATCH --mem=10G

set -euo pipefail

filenum=$SLURM_ARRAY_TASK_ID
# filenum=1

# INPUT="gLike_3G09_12273717_$filenum.log"
# INPUT="CMA_kappa_10000_12527396_$filenum.log"
# OUTPUT="../Analysis/CMA_ESES_3G09_no_m1_kappa_10000_chronological_$filenum.csv"

# python scrape_log_chronological.py "$INPUT" -o "$OUTPUT"

# echo "Done: $OUTPUT"

INPUT="CMA_true_14482288_$filenum.log"
OUTPUT="../Analysis/CMA_true_chronological_${filenum}.csv"
python scrape_log_chronological.py "$INPUT" --output "$OUTPUT" --include-initial

INPUT="CMA_true_bug_14482858_$filenum.log"
OUTPUT="../Analysis/CMA_true_bug_chronological_${filenum}.csv"
python scrape_log_chronological.py "$INPUT" --output "$OUTPUT" --include-initial

INPUT="CMA_fixed_bug_14482845_$filenum.log"
OUTPUT="../Analysis/CMA_fixed_bug_chronological_${filenum}.csv"
python scrape_log_chronological.py "$INPUT" --output "$OUTPUT" --include-initial

INPUT="CMA_fixed_14482278_$filenum.log"
OUTPUT="../Analysis/CMA_fixed_chronological_${filenum}.csv"
python scrape_log_chronological.py "$INPUT" --output "$OUTPUT" --include-initial
