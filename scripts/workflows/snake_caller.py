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

#Defining some globals
work = '/work/FAC/FBM/DMF/smitri/evomicrocomm/genome_size/data/'
s = Samples()

def submit(files):
    """Basic snakemake calling"""
    cluster_config = '--cluster-config cluster.json --cluster \
        "sbatch --mem={cluster.mem} -t {cluster.time} -c {threads}"'
    cmd = ['snakemake','--rerun-incomplete','-j','100',cluster_config,files]
    subprocess.call(' '.join(cmd),shell=True)

def snake_test():
    submit(join(work,'T22.4.5','oa','report.md'))

def strain_caller(output_file):
    """Grabbing all the directories and formatting them as snakemake wildcards."""
    output = []
    for sample in s.strains['Ochrobactrum anthropi']:
        if sample['platform'] == 'illumina':
            output.append(sample['name'])
    files = join(work,'{'+','.join(output)+'}',\
        s.abbreviations['Ochrobactrum anthropi'],output_file)
    submit(files)

if __name__ == '__main__':
    #strain_caller(join('snippy','snps.tab'))
    strain_caller('report.md')
    #snake_test()
