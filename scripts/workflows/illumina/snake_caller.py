#!/usr/bin/env python
#
# Reserve 1 CPUs for this job
#
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#
# Request it to run this for DD:HH:MM with ?G per core
#
#SBATCH --time=72:00:00
#
import subprocess
from os.path import join
from samples import Samples
import glob
import sys

"""
################################################################################
Author: https://github.com/nahanoo
This is a really usefull script to call snakemake. For snakemake you
define the output file you want to be generated. Snakemake checks all rules
and creates missing input files.
################################################################################
"""


#Defining some globals
work = '/work/FAC/FBM/DMF/smitri/evomicrocomm/genome_size/data/'
s = Samples()

def submit(files):
    """Basic snakemake calling taking the desired output file as input."""
    cluster_config = '--cluster-config cluster.json --cluster \
        "sbatch --mem={cluster.mem} -t {cluster.time} -c {threads}"'
    cmd = ['snakemake','--rerun-incomplete','-j','100',cluster_config,files]
    subprocess.call(' '.join(cmd),shell=True)

def strain_caller(strain,output_file):
    """This allows to create desired files per strain."""
    output = []
    for sample in s.strains[strain]:
        if sample['platform'] == 'illumina':
            output.append(sample['name'])
    files = join(work,'{'+','.join(output)+'}',\
        s.abbreviations[strain],output_file)
    submit(files)

def sample_caller(output_file):
    output = []
    for species,samples in s.strains.items():
        for sample in samples:
            if sample['platform'] == 'illumina':
                output.append(sample['name'])
                break
    submit(join(work,'{'+','.join(output)+'}',\
        output_file))



if __name__ == '__main__':
    """Example of a function call. Super nice that abbrevieations work.
    Example how to run this script: sbatch snake_caller.py at
    It's nice that this script then also runs as a sleeper on the cluster.
    """
    #strain_caller(s.abbreviations[sys.argv[1]],join('var.vcf'))
    for name in ['at']:
        strain_caller(s.abbreviations[name],join('spades','contigs.fasta'))
