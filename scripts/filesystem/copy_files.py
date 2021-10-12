import shutil
import glob
import pandas as pd
import os
from get_samples import Samples


class Illumina():
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

work = '/work/FAC/FBM/DMF/smitri/evomicrocomm/genome_size/data/illumina/'

def copy_reads():
    src_dir = '/nas/FAC/FBM/DMF/smitri/evomicrocomm/D2c/00_ILLUMINA_DATA/02_TRIMMOMATIC/'
    df = pd.read_csv('sample_sheet.csv')
    for sample,fname in zip(df['sample'],df['fname']):
        print(sample)
        if not os.path.exists(os.path.join(work,sample)):
            os.mkdir(os.path.join(work,sample))
        src_files = sorted(glob.glob(os.path.join(src_dir,fname+'*trimmed*')))
        shutil.copyfile(src_files[0],os.path.join(work,sample,'read1.fastq.gz'))
        shutil.copyfile(src_files[1],os.path.join(work,sample,'read2.fastq.gz'))

def copy_bams():
    src_dir = '/nas/FAC/FBM/DMF/smitri/evomicrocomm/D2c/00_ILLUMINA_DATA/04_MAPPING/04.1_BADMAPREMOVE/'

    df = pd.read_csv('sample_sheet.csv')
    for sample,fname in zip(df['sample'],df['fname']):
        files = sorted(glob.glob(os.path.join(src_dir,fname+'*.bam*')))
        shutil.copyfile(files[0],os.path.join(work,sample,'aligned.reads.sorted.bam'))
        shutil.copyfile(files[1],os.path.join(work,sample,'aligned.reads.sorted.bam.bai'))

def split_dfs():
    def create_dirs():
        for d in glob.glob(os.path.join(work,'T*')):
            if not os.path.exists(os.path.join(work,d,'at')):
                os.mkdir(os.path.join(work,d,'at'))
                os.mkdir(os.path.join(work,d,'ct'))
                os.mkdir(os.path.join(work,d,'ms'))
                os.mkdir(os.path.join(work,d,'oa'))

    strains = dict()
    strains['at'] = ['AGTU001.0001.c01','AGTU001.0001.c02','AGTU001.0001.p01',\
        'AGTU001.0001.p02','AGTU001.0001.p03']
    strains['ct'] = ['COTE001.0001.c01','COTE001.0001.c02']
    strains['ms'] = ['MISA001.0001.c01']
    strains['oa'] = ['OCAN001.0001.c01','OCAN001.0001.c02']
        
    for d in glob.glob(os.path.join(work,'T*')):
        depth = pd.read_csv(os.path.join(work,d,'depth.tsv'),sep='\t',\
            names=['chromosome','position','depth'])
        for strain in strains.keys():
            dfs = []
            for chromosome in strains[strain]:
                df = depth[depth['chromosome'] == chromosome]
                dfs.append(df)
            out = pd.concat(dfs)
            out.to_csv(os.path.join(work,d,strain,'depth.tsv'),sep='\t',index=False)
    
def copy_references():
    at = '/users/eulrich/evomicrocomm/references/at/at.fasta'
    ct = '/users/eulrich/evomicrocomm/references/ct/ct.fasta'
    ms = '/users/eulrich/evomicrocomm/references/ms/ms.fasta'
    oa = '/users/eulrich/evomicrocomm/references/oa/oa.fasta'

    for a in glob.glob(os.path.join(work,'*','at')):
        shutil.copyfile(at,a+'/reference.fasta')
    for c in glob.glob(os.path.join(work,'*','ct')):
        shutil.copyfile(ct,c+'/reference.fasta')
    for m in glob.glob(os.path.join(work,'*','ms')):
        shutil.copyfile(ms,m+'/reference.fasta')
    for o in glob.glob(os.path.join(work,'*','oa')):
        shutil.copyfile(oa,o+'/reference.fasta')