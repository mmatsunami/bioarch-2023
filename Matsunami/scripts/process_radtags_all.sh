#!/bin/bash
#SBATCH --mem=64G 

SAMPLE=$1

module load singularity/3.8.3

singularity exec -B /lustre8,/home /usr/local/biotools/s/stacks:2.65--hdcf5f25_0 \
process_radtags -P \
-1 /home/bioarchaeology-pg/data/11/rawdata/${SAMPLE}_R1.fq.gz \
-2 /home/bioarchaeology-pg/data/11/rawdata/${SAMPLE}_R2.fq.gz \
-o samples \
-c -q --renz_1 ecoRI --renz_2 mseI

mv samples/${SAMPLE}_R1.1.fq.gz samples/${SAMPLE}.1.fq.gz
mv samples/${SAMPLE}_R2.2.fq.gz samples/${SAMPLE}.2.fq.gz
