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

metaphlan --nproc 32 $work/$1/read1.fq.gz,$work/$1/read2.fq.gz --input_type fastq \
-s $work/$1/$1.sam.bz2 --bowtie2out $work/$1/$1.bowtie2.bz2 \
-o $work/$1/$1._profiled.tsv

sample2markers.py -i $work/$1/$1.sam.bz2 -o $work/$1/$1.pkl -n 32
