import pandas as pd
from os.path import join
from samples import Samples
from Bio import SeqIO
import plotly.express as px
import numpy as np
s = Samples()


"""This is a short script that annotates variants for Illumina data."""

def annotate_variants(abb):
    """Prints annotations for identified SNPs."""
    def annotate_pos(gbk, c, p):
        for (start, end), (gene, product) in gbk[c].items():
            if p in range(start, end):
                return [(start, end), gene, product]
        return False

    f = join('..', 'annotations', abb + '_renamed.tsv')
    df = pd.read_csv(f, sep='\t')
    # For plotting we renamed contigs to at_0 etc.
    # Rename of contigs in annotations for hashing.
    contigs = {c: abb+'_'+str(i)
               for i, c in enumerate(sorted(set(df['Sequence Id'])))}
    for i, chrom in enumerate(df['Sequence Id']):
        df.at[i, 'Sequence Id'] = contigs[chrom]
    gbk = {contig: {} for contig in sorted(set(df['Sequence Id']))}
    for i, row in df.iterrows():
        if pd.isna(row['Gene']):
            gene = 'Unknown'
        else:
            gene = row['Gene']
        if pd.isna(row['Product']):
            product = 'hypothetical protein'
        else:
            product = row['Product']
        gbk[row['Sequence Id']][(row['Start'], row['Stop'])] = (gene, product)

    in_files = [join('..', 'variants', 'variants_comp_mapping.csv'), join(
        '..', 'variants', 'snps_freebayes_comp_mapping.csv'),
        join('..', 'variants', 'snps_pacbio_renamed.csv')]
    out_files = [join('..', 'annotations', abb + '_variants_annotations.csv'),
                 join('..', 'annotations', abb + '_snps_annotations.csv'),
                 join('..', 'annotations', abb + '_pacbio_annotations.csv')]
    for j, in_file in enumerate(in_files):
        df = pd.read_csv(in_file)
        df = df[df['strain'] == s.abbreviations[abb]]
        df.insert(len(df.columns), 'gene_start', None)
        df.insert(len(df.columns), 'gene_end', None)
        df.insert(len(df.columns), 'gene', None)
        df.insert(len(df.columns), 'product', None)
        df.insert(len(df.columns), 'sequence', None)
        for i, row in df.iterrows():
            a = annotate_pos(gbk, row['chrom'], row['pos'])
            if a:
                pass
            else:
                a = [('Not annotated', 'Not annotated'),
                     'Not annotated', 'Not annotated']
            df.at[i, 'gene_start'], df.at[i, 'gene_end'] = a[0][0], a[0][1]
            df.at[i, 'gene'], df.at[i, 'product'] = a[1], a[2]
        contigs = {contig.name: str(contig.seq) for contig in SeqIO.parse(
            join('..', 'annotations', abb+'.fasta'), 'fasta')}
        for i, row in df.iterrows():
            if row['gene'] != 'Not annotated':
                seq = contigs[row['chrom']][int(
                    row['gene_start']):int(row['gene_end'])]
            else:
                seq = 'Not annotated'
            df.at[i, 'sequence'] = seq
        df.to_csv(out_files[j], index=False)

