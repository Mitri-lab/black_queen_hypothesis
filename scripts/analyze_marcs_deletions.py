from samples import Samples
import os
import pandas as pd
from Bio import SeqIO
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import json

"""
################################################################################
Author: https://github.com/nahanoo
This script is not important anymore. It analyses deletions identified by a previous
analysis. We wanted to check, if the identified deletions have a reasonable sequence
composition. Unfortunately we saw, that most deletions have extremely high,
or extremely low GC conten. This is visualized by plotting the GC conten of a
deleted sequence and comparing it to the histogram showing the GC content of every
150-mer.
################################################################################
"""

#Globally initiating samples class
s = Samples()

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
    work = '/work/FAC/FBM/DMF/smitri/evomicrocomm/genome_size/data/'
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
    with open('gc_contents.json','r') as handle:
        data = handle.read()
    gc_contents = json.loads(data)

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