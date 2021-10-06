from get_samples import Samples
import os
import pandas as pd
from Bio import SeqIO
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import json

s = Samples()

def get_gc_content(sequence):
    return 100.0*len([base for base in sequence if base in "GC"])/len(sequence)

def get_gc_contents():
    window_size = 150
    gc_contents = dict()
    for strain in s.strains:
        print(strain)
        gc_content = []
        reference = {contig.name:contig for contig in SeqIO.parse(s.references[strain],'fasta')}
        for chromosome,record in reference.items():
            for position in range(len(record)-window_size):
                gc_content.append(get_gc_content(record[position:position+window_size]))
        gc_contents[strain] = gc_content
    return gc_contents

def gc_histogram(fig,gc_contents,strain):
    fig.append_trace(go.Histogram(x=gc_contents[strain],name='gc content',nbinsx=20),1,1)
    fig.update_layout(
        xaxis_title = 'GC content in %',
        yaxis_title = 'counts'
    )
    return fig

def deletion_histogram(fig,strain):
    work = '/work/FAC/FBM/DMF/smitri/evomicrocomm/genome_size/data/'
    gc_contents = []
    for sample in s.strains[strain]:
        f = os.path.join(work,sample['dir'],'grouped_deletions.tsv')
        contigs = {contig.name:contig for contig in SeqIO.parse(s.references[strain],'fasta')}
        if os.path.exists(f):
            df = pd.read_csv(f,sep='\t')
            for chromosome,position,length in zip(df['chromosome'],df['position'],df['length']):
                sequence = str(contigs[chromosome][position-50:position+length+50].seq)
                if len(sequence) > 0:
                    gc_contents.append(get_gc_content(sequence))
    fig.append_trace(go.Histogram(x=gc_contents,name='gc content',nbinsx=20),1,2)
    return fig


with open('gc_contents.json','r') as handle:
    data = handle.read()
gc_contents = json.loads(data)


for strain in s.strains:
    print(strain)
    fig = make_subplots(1,2,subplot_titles=("Reference","Deletions"))
    fig = gc_histogram(fig,gc_contents,strain)
    fig = deletion_histogram(fig,strain)
    fig.update_layout(
        title = strain,
        xaxis_title = 'GC content in %',
        yaxis_title = 'Counts'
    )
    fig.update_traces(showlegend=False)
    fig.write_image(strain.replace(' ','_')+'.png')
