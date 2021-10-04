#!/bin/sh
# Reserve 24 CPUs for this job
#
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#
# Request it to run this for HH:MM:SS with ?G per core
#
#SBATCH --time=00:30:00
#
minimap2 -t 16 --secondary=no -ax asm20 $1/reference.fasta $1/corrected_reads.fasta > $1/alignment.minimap.sam
samtools view -S -b $1/alignment.minimap.sam > $1/alignment.minimap.bam
samtools sort $1/alignment.minimap.bam -o $1/alignment.minimap.sorted.bam
samtools index $1/alignment.minimap.sorted.bam
samtools depth -aa $1/alignment.minimap.sorted.bam -o $1/depth.tsv
detect_deletions --output_no_alignment_regions $1/alignment.minimap.sorted.bam $1/
