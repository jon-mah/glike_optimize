#!/bin/bash
#SBATCH --array=1-10

set -euo pipefail

filenum=$SLURM_ARRAY_TASK_ID

INPUT="gLike_3G09_12273717_$filenum.log"
OUTPUT="../Analysis/gLike_3G09_chronological_$filenum.csv"

python scrape_log_chronological.py "$INPUT" -o "$OUTPUT"

echo "Done: $OUTPUT"
