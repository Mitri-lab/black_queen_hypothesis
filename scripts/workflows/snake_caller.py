import glob
import subprocess
from os.path import join

#Defining some globals
work = '/work/FAC/FBM/DMF/smitri/evomicrocomm/genome_size/data/'
cluster_config = '--cluster-config cluster.json --cluster "sbatch --mem={cluster.mem} -t {cluster.time} -c {threads}"'

def submit(output):
    cmd = ['snakemake','-p','-j','100',cluster_config,output]
    subprocess.call(' '.join(cmd),shell=True)

def get_at_dirs():
    

output = join(work,'T33.3.5','ct','mapped_reads.sam')
cmd = ['snakemake','-p','-j','100',cluster_config,output]
subprocess.call(' '.join(cmd),shell=True)