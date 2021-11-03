from samples import Samples
from plotting import Plotting
from analyze_marcs_deletions import get_gc_content
import pandas as pd
from os.path import join
from os.path import exists
import math
from Bio import SeqIO
import glob

"""
################################################################################
Author: https://github.com/nahanoo
This script analyzes all outputs generated from illumina sample processing.
If you are interested in the sample processing you can check out the Snakemake
workflow: https://github.com/nahanoo/black_queen_hypothesis/blob/main/scripts/workflows/illumina/Snakefile.
All plots call this small plotting class in plotting.py
################################################################################
"""

s = Samples()
p = Plotting()

def plot_snps():
    """Plots n SNPs of every strain per treatment."""
    for strain in s.strains:
        #Get all treatments of a strain
        treatments = s.treatments[strain]
        #df for storing n_snps
        n_snps = pd.DataFrame(columns=treatments,index=[sample['name'] \
            for sample in s.strains[strain] if sample['platform']== 'illumina'])
        #Iterating over every sample of a strain
        for sample in s.strains[strain]:
            if sample['platform'] == 'illumina':
                snps = join(sample['dir_name'],'snippy','snps.tab')
                if exists(snps):
                    #Storing found n snps in df
                    n_snps.at[sample['name'],sample['treatment']] = \
                        len(pd.read_csv(snps,sep='\t'))

        #Plotting n snps
        fig = p.subplot_treatments(strain,n_snps)
        #Updating plot labels and dumping to png
        title = 'N SNPs in '+strain
        fig.update_layout(
            xaxis_title='samples',
            yaxis_title='n SNPs',
            width=len(n_snps) * 30,
            title=title)
        fig.update_traces(showlegend=False)
        fig.write_image(join('..','plots','snps',title.replace(' ','_')+'.png'))

def plot_genes():
    """This plots the counts of a mutated gene over the time series
    of an experiment. This visualizes in how many microcosms we see
    a mutation arise."""
    #Iterating over strains
    for strain,samples in s.strains.items():
        treatments = s.treatments[strain]
        samples = [sample for sample in samples if sample['platform'] == 'illumina']
        #We create a plot per treatment, which is why we iterate over all treatments.
        for treatment in treatments:
            timepoints = ['T11','T22','T33','T44']
            out = pd.DataFrame(columns=timepoints)
            for sample in samples:
                #We only want matching treatments
                if sample['treatment'] == treatment:
                    snps = join(sample['dir_name'],'snippy','snps.tab')
                    if exists(snps):
                        for gene in set(pd.read_csv(snps,sep='\t').dropna()['GENE']):
                            #If gene not in index we neet to set value to 1 else we
                            #add 1 to the count
                            if gene in out['treatment'].dropna().index:
                                out.at[gene,sample['timepoint']] += 1
                            else:
                                out.at[gene,sample['timepoint']] = 1

            #Plotting gene counts
            fig = p.trajectories(out)
            #Changing labels and dumping to png
            title = ['Mutated','genes','in',strain,'in','treatment',str(treatment)]
            fig.update_layout(
                    title = ' '.join(title),
                    xaxis_title = 'timepoints',
                    yaxis_title = 'observed in n microcosms'
                )
            fig.write_image(join('..','plots','genes',' '.join(title).replace(' ','_')+'.png'))

def plot_products():
    """Plotting in how many samples a product was mutated."""
    #Iterating over all strains
    for strain,samples in s.strains.items():
        treatments = s.treatments[strain]
        #Setting up df
        n_products = pd.DataFrame(columns=treatments)
        for sample in samples:
            if sample['platform'] == 'illumina':
                f = join(sample['dir_name'],'snippy','snps.tab')
                if exists(f):
                    products = set(pd.read_csv(f,sep='\t').dropna()['PRODUCT'])
                    for product in products:
                        #Checking if product is in df already, if so adding 1
                        if product in n_products[sample['treatment']].dropna().index:
                            n_products.at[product,sample['treatment']] += 1
                        else:
                            n_products.at[product,sample['treatment']] = 1

        #Creating plo
        fig = p.subplot_products(strain,n_products)
        title = ['Products','affected','by','mutations','in',strain]
        #Updating labels and dumping plot to png
        fig.update_layout(overwrite=True,
                title = ' '.join(title),
                height = len(n_products) * 50
            )
        fig.update_xaxes(title_text='observed mutated product n times',row=len(treatments),col=1)
        fig.update_yaxes(title_text='products',row=int(math.ceil(len(treatments)/2)),col=1)
        fig.update_traces(showlegend=False)
        fig.write_image(join('..','plots','products',' '.join(title).replace(' ','_')+'.png'))

def write_gc_content():
    """This is a little helper function to get the GC content of the wild-type genomes"""
    #Setting up df and iterating over all strains
    df = pd.DataFrame(columns=['gc content'])
    for strain in s.strains:
        sequence = str()
        reference = s.references[strain]
        contigs = [contig for contig in SeqIO.parse(reference,'fasta')]
        for contig in contigs:
            #Adding sequence as string
            sequence += contig.seq
        #Getting gc content of sequence
        gc_content = get_gc_content(sequence)
        #Writing to df and to file
        df.at[strain,'gc content'] = gc_content
    df.index.name = 'strain'
    fname = join('..','tables','gc_content','gc_content.csv')
    df.to_csv(fname)