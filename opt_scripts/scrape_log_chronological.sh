#!/bin/bash
#SBATCH --array=1-10

set -euo pipefail

filenum=$SLURM_ARRAY_TASK_ID

# INPUT="gLike_3G09_12273717_$filenum.log"
INPUT="CMA_kappa_10000_12527396_$filename.log"
OUTPUT="../Analysis/CMA_ESES_3G09_no_m1_kappa_10000_chronological_$filenum.csv"

python scrape_log_chronological.py "$INPUT" -o "$OUTPUT"

echo "Done: $OUTPUT"
