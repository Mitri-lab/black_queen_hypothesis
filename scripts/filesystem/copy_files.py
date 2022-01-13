import shutil
import glob
import pandas as pd
import os
from os.path import join,exists
from samples import Samples

def create_dirs():
    """This creates the directory sturcture for the hgt analysis of the pacbio data"""
    s = Samples()
    for ancestral, samples in s.strains.items():
        for sample in samples:
            if sample["platform"] == "pacbio":
                p = join(sample["dir_name"], "hgt")
                if not exists(p):
                    mkdir(p)
                for strain in s.strains_per_treatment[sample["treatment"]]:
                    reference = join(
                        sample["dir_name"], "hgt", s.abbreviations[strain] + ".fasta"
                    )
                    if not exists(reference):
                        symlink(
                            s.references[strain],
                            reference,
                        )


class Illumina():
    """This is quick code which allows to restore all files from NAS.
    Those moved files can then be used for the Snakemkae workflow
    https://github.com/nahanoo/black_queen_hypothesis/tree/main/scripts/workflows/illumina
    """
    def __init__(self):
        self.target_dir = '/work/FAC/FBM/DMF/smitri/evomicrocomm/genome_size/data/'
        self.src_dir = '/nas/FAC/FBM/DMF/smitri/evomicrocomm/D2c/00_ILLUMINA_DATA/02_TRIMMOMATIC/'
        self.samples = Samples()

    def create_dirs(self):
        """Iterating over all Illumina samples and creating directories"""
        for strain in self.samples.strains:
            for sample in self.samples.strains[strain]:
                if sample['platform'] == 'illumina':
                    if not os.path.exists(sample['dir']):
                        os.makedirs(sample['dir'])
    
    def copy_references(self):
        """Copying references to every directory. Bit overkill but convenient."""
        for strain in self.samples.strains:
            for sample in self.samples.strains[strain]:
                if sample['platform'] == 'illumina':
                    if not os.path.exists(os.path.join(sample['dir'],'reference.fasta')):
                        shutil.copyfile(self.samples.references[strain],os.path.join(sample['dir'],'reference.fasta'))

    def copy_reads(self):
        df = pd.read_csv('illumina_names.csv')
        for fname,sample in zip(df['fname'],df['sample']):
            reads = sorted(glob.glob(os.path.join(self.src_dir,fname+'*trimmed.fq.gz')))
            names = ['read1.fq.gz','read2.fq.gz']
            for read,name in zip(reads,names):
                cmd = ['sbatch','copy_files.sh',read,\
                    os.path.join(self.target_dir,sample,name)]
                if not os.path.exists(os.path.join(self.target_dir,sample,name)):
                    print(sample)
                    shutil.copyfile(read,os.path.join(self.target_dir,sample,name))