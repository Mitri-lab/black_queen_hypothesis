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

work = '/work/FAC/FBM/DMF/smitri/evomicrocomm/genome_size/data/ancestors'
s = Samples()


def submit(files):
    """Basic snakemake calling taking the desired output file as input."""
    cluster_config = '--cluster-config cluster.json --cluster \
        "sbatch --mem={cluster.mem} -t {cluster.time} -c {threads}"'
    cmd = ['snakemake', '--rerun-incomplete','--cores 16', files]
    print(cmd)
    subprocess.call(' '.join(cmd), shell=True)


def caller(output_file):
    output = []
    for strain in s.strains:
        output.append(s.abbreviations[strain])
    files = join(work, '{'+','.join(output)+'}', output_file)
    submit(files)


caller(join('bakta','assembly.contigs.polypolish.gbff'))
