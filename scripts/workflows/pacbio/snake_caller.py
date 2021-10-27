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

#Defining some globals
work = '/work/FAC/FBM/DMF/smitri/evomicrocomm/genome_size/data/'
s = Samples()

def submit(files):
    """Basic snakemake calling"""
    cluster_config = '--cluster-config cluster.json --cluster \
        "sbatch --mem={cluster.mem} -t {cluster.time} -c {threads}"'
    cmd = ['snakemake','--rerun-incomplete','-j','100',cluster_config,files]
    subprocess.call(' '.join(cmd),shell=True)

def single_caller(path):
    submit(path)

def all_caller(output_file):
    output = []
    for strain,samples in s.strains.items():
        for sample in samples:
            if sample['platform'] == 'pacbio':
                output.append(sample['name'])
    files = join(work,'{'+','.join(output)+'}',output_file)
    submit(files)

if __name__ == '__main__':
    #single_caller(sys.argv[1])
    all_caller(join('prokka','prokka.gbk'))