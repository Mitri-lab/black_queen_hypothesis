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
plot_gc_bias --plotting_only $1/reference.fasta $1/alignment.minimap.sorted.bam $1/labels.json $1/