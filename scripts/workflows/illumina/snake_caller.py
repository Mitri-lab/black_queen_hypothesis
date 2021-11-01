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

def strain_caller(strain,output_file):
    """Grabbing all the directories and formatting them as snakemake wildcards."""
    output = []
    for sample in s.strains[strain]:
        if sample['platform'] == 'illumina':
            output.append(sample['name'])
    files = join(work,'{'+','.join(output)+'}',\
        s.abbreviations[strain],output_file)
    submit(files)

def time_point_caller(output_file):
    dirs = [d.split('/')[-1] for d in glob.glob(join(work,'T44.4.*'))]
    submit(join(work,'{'+','.join(dirs)+'}','{at,ct,oa,ms}',output_file))

if __name__ == '__main__':
    #strain_caller(sys.argv[1],join('snippy','snps.tab'))
    strain_caller(s.abbreviations[sys.argv[1]],join('snippy','snps.tab'))
    #snake_test()
    #time_point_caller('report.md')
    #single_caller(sys.argv[1])
    #strain_caller(s.abbreviations[sys.argv[1]],'flagstat.tsv')
