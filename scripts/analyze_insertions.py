from Bio import SeqIO
import pandas as pd
from glob import glob
from os.path import join


def parse_insertions(tsv):
    return pd.read_csv(
        tsv, sep="\t", usecols=["chromosome", "position", "length"]
    ).drop_duplicates()


def get_contigs(fasta):
    contigs = [contig for contig in SeqIO.parse(fasta, "fasta")]
    c = {contig.id: len(contig) for contig in contigs}
    return c


def create_df(fasta, insertions):
    contigs = get_contigs(fasta)
    insertions.insert(3, "contig_length", None)
    for i, row in insertions.iterrows():
        insertions.at[i, "contig_length"] = contigs[row.chromosome]
    return insertions


def all():
    work = '/users/eulrich/work/genome_size/data'
    tsvs = glob(join(work,'*','mutant_to_parent.noalignments.tsv'))
    for tsv in tsvs:
        fasta = tsv.replace('mutant_to_parent.noalignments.tsv','assembly.fasta')
        insertions = parse_insertions(tsv)
        df = create_df(fasta, insertions)
        print(tsv.split('/')[-2],'##########################')
        print(df)

def single():
    work = '/users/eulrich/work/genome_size/data'
    sample = 'Ct43.1'
    tsv = join(work,sample,'mutant_to_parent.noalignments.tsv')
    fasta = join(work,sample,'assembly.fasta')
    insertions = parse_insertions(tsv)
    df = create_df(fasta, insertions)
    print(df)