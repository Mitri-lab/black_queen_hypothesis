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
from subprocess import call
from os.path import join
from os.path import exists
from samples import Samples
from os import mkdir

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
    """Basic snakemake calling"""
    cluster_config = '--cluster-config cluster.json --cluster \
        "sbatch --mem={cluster.mem} -t {cluster.time} -c {threads}"'
    cmd = ['snakemake','--rerun-incomplete','-j','100',cluster_config,files]
    call(' '.join(cmd),shell=True)

def metaphlan():
    """Basic snakemake calling taking the desired output file as input."""
    outputs = []
    for strain,samples in s.strains.items():
        for sample in samples:
            if sample['platform'] == 'illumina':
                output = join(sample['name'],sample['name']+'.pkl')
                outputs.append(output)
    files = join(work,'{'+','.join(outputs)+'}')
    submit(files)

def db_markers():
    for strain in s.strains.keys():
        call(['extract_markers.py','-c','s__'+strain.replace(' ','_')\
            ,'-o',work+s.abbreviations[strain]])

if __name__ == '__main__':
    """Example how to run this script: sbatch snake_caller.py at
    It's nice that this script then also runs as a sleeper on the cluster.
    """
    #metaphlan()
    pass