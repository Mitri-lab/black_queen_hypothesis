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
    subprocess.call(' '.join(cmd),shell=True)


def all_caller(output_file):
    """Basic snakemake calling taking the desired output file as input."""
    output = []
    for strain,samples in s.strains.items():
        for sample in samples:
            if sample['platform'] == 'pacbio':
                output.append(sample['name'])
    files = join(work,'{'+','.join(output)+'}',output_file)
    submit(files)

if __name__ == '__main__':
    """Example how to run this script: sbatch snake_caller.py at
    It's nice that this script then also runs as a sleeper on the cluster.
    """
    #submit(join(work,'At42.1','deletions.annotated.tsv'))
    all_caller(join('deletions.annotated.tsv'))