from samples import Samples
from plotting import Plotting
from analyze_marcs_deletions import get_gc_content
import pandas as pd
from os.path import join
from os.path import exists
import math
from Bio import SeqIO
import glob

s = Samples()
p = Plotting()

def plot_snps():
    for strain in s.strains:
        #Get all treatments of a strain
        treatments = s.treatments[strain]
        #Create subplot titles
        n_snps = pd.DataFrame(columns=treatments,index=[sample['name'] \
            for sample in s.strains[strain] if sample['platform']== 'illumina'])
        for counter,treatment in enumerate(treatments):
            #Get all sample names of a treatment
            for sample in s.strains[strain]:
                #Getting genome lenght and n contigs
                if sample['platform'] == 'illumina':
                    snps = join(sample['dir_name'],'snippy','snps.tab')
                    if exists(snps):
                        n_snps.at[sample['name'],sample['treatment']] = \
                            len(pd.read_csv(snps,sep='\t'))
        fig = p.subplot_treatments(strain,n_snps)
        title = 'N SNPs in '+strain
        fig.update_layout(
            xaxis_title='samples',
            yaxis_title='n SNPs',
            width=len(n_snps) * 30,
            title=title)
        fig.update_traces(showlegend=False)
        fig.write_image(join('..','plots','snps',title.replace(' ','_')+'.png'))

def plot_genes():
    """This plots the mutated genes in n microocms."""
    for strain,samples in s.strains.items():
        treatments = s.treatments[strain]
        samples = [sample for sample in samples if sample['platform'] == 'illumina']
        for treatment in treatments:
            genes = []
            for sample in samples:
                if sample['treatment'] == treatment:
                    snps = join(sample['dir_name'],'snippy','snps.tab')
                    if exists(snps):
                        genes += pd.read_csv(snps,sep='\t').dropna()['GENE'].to_list()
            timepoints = ['T11','T22','T33','T44']
            out = pd.DataFrame(0,columns=timepoints,index=set(genes))
            for sample in samples:
                if sample['treatment'] == treatment:
                    snps = join(sample['dir_name'],'snippy','snps.tab')
                    if exists(snps):
                        for gene in set(pd.read_csv(snps,sep='\t').dropna()['GENE']):
                            out.at[gene,sample['timepoint']] += 1 
            fig = p.trajectories(out)
            title = ['Mutated','genes','in',strain,'in','treatment',str(treatment)]
            fig.update_layout(
                    title = ' '.join(title),
                    xaxis_title = 'timepoints',
                    yaxis_title = 'observed in n microcosms'
                )
            fig.write_image(join('..','plots','genes',' '.join(title).replace(' ','_')+'.png'))

def plot_products():
    for strain,samples in s.strains.items():
        treatments = s.treatments[strain]
        n_products = pd.DataFrame(columns=treatments)
        for sample in samples:
            if sample['platform'] == 'illumina':
                f = join(sample['dir_name'],'snippy','snps.tab')
                if exists(f):
                    products = set(pd.read_csv(f,sep='\t').dropna()['PRODUCT'])
                    for product in products:
                        if product in n_products[sample['treatment']].dropna().index:
                            n_products.at[product,sample['treatment']] += 1
                        else:
                            n_products.at[product,sample['treatment']] = 1
        fig = p.subplot_products(strain,n_products)
        title = ['Products','affected','by','mutations','in',strain]
        fig.update_layout(overwrite=True,
                title = ' '.join(title),
                height = len(n_products) * 50
            )
        fig.update_xaxes(title_text='observed mutated product n times',row=len(treatments),col=1)
        fig.update_yaxes(title_text='products',row=int(math.ceil(len(treatments)/2)),col=1)
        fig.update_traces(showlegend=False)
        fig.write_image(join('..','plots','products',' '.join(title).replace(' ','_')+'.png'))

def write_gc_content():
    df = pd.DataFrame(columns=['gc content'])
    for strain in s.strains:
        sequence = str()
        reference = s.references[strain]
        contigs = [contig for contig in SeqIO.parse(reference,'fasta')]
        for contig in contigs:
            sequence += contig.seq
        gc_content = get_gc_content(sequence)
        df.at[strain,'gc content'] = gc_content
    df.index.name = 'strain'
    fname = join('..','tables','gc_content','gc_content.csv')
    df.to_csv(fname)
