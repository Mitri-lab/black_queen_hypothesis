#!/bin/sh
# Reserve 24 CPUs for this job
#
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#
# Request it to run this for HH:MM:SS with ?G per core
#
#SBATCH --time=03:00:00
#
work=/work/FAC/FBM/DMF/smitri/evomicrocomm/genome_size/data
strainphlan -s $work/$1/$1.pkl -m $work/$2/s__$3.fna -r $work/$2/reference.fasta \
 -o $work/$2/strainphlan -n 8 -c s__$2 --mutation_rates

