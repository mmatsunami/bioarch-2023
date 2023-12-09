#!/bin/bash
#SBATCH --mem=64G 

SAMPLE=$1

module load singularity/3.8.3

#mapping
singularity exec -B /lustre8,/home /usr/local/biotools/b/bwa:0.7.8--hed695b0_5 \
bwa mem /home/bioarchaeology-pg/data/11/reference/yaponesia_reference.fasta \
samples/${SAMPLE}.1.fq.gz \
samples/${SAMPLE}.2.fq.gz \
> mapping/${SAMPLE}.sam

#sam2bam
singularity exec -B /lustre8,/home /usr/local/biotools/s/samtools:1.18--hd87286a_0 \
samtools view -bS mapping/${SAMPLE}.sam > mapping/${SAMPLE}.bam

#sort
singularity exec -B /lustre8,/home /usr/local/biotools/s/samtools:1.18--hd87286a_0 \
samtools sort mapping/${SAMPLE}.bam -o mapping/${SAMPLE}.bam
