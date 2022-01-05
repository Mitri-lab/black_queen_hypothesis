from Bio import SeqIO
import pandas as pd
from glob import glob
from os.path import join


def get_contigs(fasta):
    contigs = [contig for contig in SeqIO.parse(fasta, "fasta")]
    c = {contig.id: contig.seq for contig in contigs}
    return c


def get_sequences():
    work = '/users/eulrich/work/genome_size/data'
    df = pd.read_csv("sequences_of_interest.csv")
    df.insert(5,'sequence',None)
    for i, row in df.iterrows():
        c = row.chromosome
        p = row.position
        l = row.length
        s = row.sample_name
        fasta = join(work,s,'assembly.fasta')
        contigs = get_contigs(fasta)
        df.at[i,'sequence'] = str(contigs[c][p:p+l])
    return df