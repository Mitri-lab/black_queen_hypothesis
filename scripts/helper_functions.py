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

"""Functions to analyze GC content of deletions"""
def get_gc_content(sequence):
    """This function returns gc content for string sequence"""
    return 100.0*len([base for base in sequence if base in "GC"])/len(sequence)

def get_gc_contents():
    """This function stores the GC contents of every 150-mer of the reference sequences."""
    window_size = 150
    gc_contents = dict()
    for strain in s.strains:
        gc_content = []
        #Parsing reference sequence by making use of Samples class
        reference = {contig.name:contig for contig in SeqIO.parse(s.references[strain],'fasta')}
        #Iterating over every contig
        for chromosome,record in reference.items():
            #Iterating over every position of the contig and caclulating gc content for window.
            for position in range(len(record)-window_size):
                gc_content.append(get_gc_content(record[position:position+window_size]))
        gc_contents[strain] = gc_content
    return gc_contents

def gc_histogram(fig,gc_contents,strain):
    """Takes a plotly subplot as input and plots gc content as histogram."""
    fig.append_trace(go.Histogram(x=gc_contents[strain],name='gc content',nbinsx=20),1,1)
    fig.update_layout(
        xaxis_title = 'GC content in %',
        yaxis_title = 'counts'
    )
    return fig

def deletion_histogram(fig,strain):
    """Parses all deletions identified in Illumina. Calculates the GC
    content of the location of the deletion. Since many deletions are very small
    50 basepairs are added to both sides."""
    #Data path
    gc_contents = []
    for sample in s.strains[strain]:
        #Not all samples have deletions thus the os.path.exists check
        f = os.path.join(sample['dir'],'grouped_deletions.tsv')
        contigs = {contig.name:contig for contig in SeqIO.parse(s.references[strain],'fasta')}
        if os.path.exists(f):
            df = pd.read_csv(f,sep='\t')
            for chromosome,position,length in zip(df['chromosome'],df['position'],df['length']):
                sequence = str(contigs[chromosome][position-50:position+length+50].seq)
                #Iterating over every deletions
                if len(sequence) > 0:
                    gc_contents.append(get_gc_content(sequence))
    #Creating histogram
    fig.append_trace(go.Histogram(x=gc_contents,name='gc content',nbinsx=20),1,2)
    return fig

def main():
    """Parses deletions, gets reference gc contents and creates subplots."""
    #I dropped gc contents of the references as json to save time.
    gc_contents = get_gc_contents()

    #Creating figures for every strain
    for strain in s.strains:
        print(strain)
        #Creating subplots
        fig = make_subplots(1,2,subplot_titles=("Reference","Deletions"))
        #Adding reference gc content trace
        fig = gc_histogram(fig,gc_contents,strain)
        #Adding deletions gc content trace
        fig = deletion_histogram(fig,strain)
        fig.update_layout(
            title = strain,
            xaxis_title = 'GC content in %',
            yaxis_title = 'Counts'
        )
        fig.update_traces(showlegend=False)
        #Write to plots directory
        fig.write_image('../plots/gc_content_deletions/'+strain.replace(' ','_')+'.png')