#!/bin/bash
#SBATCH --mem=64G 

SAMPLE=$1

module load singularity/3.8.3

singularity exec -B /lustre8,/home /usr/local/biotools/s/stacks:2.65--hdcf5f25_0 \
ustacks -t gzfastq \
-f samples/${SAMPLE}.1.fq.gz \
-o denovo_map \
-i 1 --name ${SAMPLE} -M 5 -m 3
