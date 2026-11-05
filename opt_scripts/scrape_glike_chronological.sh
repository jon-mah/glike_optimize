#!/bin/bash
#SBATCH --array=1-10

set -euo pipefail

filenum=$SLURM_ARRAY_TASK_ID
# filenum=1

# INPUT="gLike_3G09_12273717_$filenum.log"
# INPUT="CMA_kappa_10000_12527396_$filenum.log"
# OUTPUT="../Analysis/CMA_ESES_3G09_no_m1_kappa_10000_chronological_$filenum.csv"

# python scrape_log_chronological.py "$INPUT" -o "$OUTPUT"

# echo "Done: $OUTPUT"

INPUT="gL_kappa_10000_13465827_$filenum.log"
OUTPUT="../Analysis/glike_10000_chronological_${filenum}.csv"
python scrape_glike_chronological.py "$INPUT" --output "$OUTPUT"

INPUT="gL_kappa_20000_13465828_$filenum.log"
OUTPUT="../Analysis/glike_20000_chronological_${filenum}.csv"
python scrape_glike_chronological.py "$INPUT" --output "$OUTPUT"
