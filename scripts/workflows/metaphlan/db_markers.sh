#!/bin/sh
# Reserve 24 CPUs for this job
#
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#
# Request it to run this for HH:MM:SS with ?G per core
#
#SBATCH --time=00:30:00
#
work=/work/FAC/FBM/DMF/smitri/evomicrocomm/genome_size/data
extract_markers.py -c s__$1 -o $work/$2/

