from samples import Samples
import pandas as pd
from os.path import join,exists
from Bio import SeqIO
import re



"""
################################################################################
Author: https://github.com/nahanoo
This is a collection of functions for illumina data that are not relevant for the paper
but worth keeping.
################################################################################
"""

s = Samples()


def plot_effects():
    """This plots the summed effects grouped per treatment per species
    """
    snps = {strain: None for strain in s.strains}
    for strain, samples in s.strains.items():
        effects = []
        for sample in samples:
            if (sample['platform'] == 'illumina') & (sample['timepoint'] == 'T44'):
                f = join(sample['dir_name'], 'snippy', 'snps.tab')
                df = pd.read_csv(f, sep='\t').drop_duplicates()
                df = df[df['EFFECT'].notna()]
                for effect in df['EFFECT']:
                    effects.append(effect.split(' ')[0])
        effects = list(set(effects))
        columns = s.treatments[strain]
        snp = pd.DataFrame(columns=columns, index=effects)
        for sample in samples:
            if (sample['platform'] == 'illumina') & (sample['timepoint'] == 'T44'):
                f = join(sample['dir_name'], 'snippy', 'snps.tab')
                df = pd.read_csv(f, sep='\t').drop_duplicates()
                df = df[df['EFFECT'].notna()]
                effects = []
                for effect in df['EFFECT']:
                    effects.append(effect.split(' ')[0])

                for effect in set(effects):
                    mask = []
                    for e in df['EFFECT']:
                        if re.search(effect, e, flags=re.IGNORECASE):
                            mask.append(True)
                        else:
                            mask.append(False)
                    if pd.isna(snp.at[effect, sample['treatment']]):
                        snp.at[effect, sample['treatment']] = len(df[mask])
                    else:
                        snp.at[effect, sample['treatment']] += len(df[mask])
        snps[strain] = snp
    return snps


def write_gc_content():
    """Small function to get gc content of references."""
    def get_gc_content(sequence):
        """This function returns gc content for string sequence"""
        return 100.0*len([base for base in sequence if base in "GC"])/len(sequence)
    """This is a little helper function to get the GC content of the wild-type genomes"""
    # Setting up df and iterating over all strains
    df = pd.DataFrame(columns=['gc content'])
    for strain in s.strains:
        sequence = str()
        reference = s.references[strain]
        contigs = [contig for contig in SeqIO.parse(reference, 'fasta')]
        for contig in contigs:
            # Adding sequence as string
            sequence += contig.seq
        # Getting gc content of sequence
        gc_content = get_gc_content(sequence)
        # Writing to df and to file
        df.at[strain, 'gc content'] = gc_content
    df.index.name = 'strain'
    fname = join('..', 'tables', 'gc_content', 'gc_content.csv')
    df.to_csv(fname)
