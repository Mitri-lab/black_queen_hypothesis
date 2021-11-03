from samples import Samples
from samples import Experiment
from plotting import Plotting
from os.path import join
from os.path import exists
import pandas as pd
from Bio import SeqIO

"""
################################################################################
Author: https://github.com/nahanoo
This script analyzes all outputs generated from PacBio sample processing.
If you are interested in the sample processing you can check out the Snakemake
workflow: https://github.com/nahanoo/black_queen_hypothesis/blob/main/scripts/workflows/pacbio/Snakefile.
All plots call this small plotting class in plotting.py
################################################################################
"""

s = Samples()
p = Plotting()
e = Experiment()

def plot_genome_length():
    """This function plots the assembly length and the number of contigs
    produced by assembling the PacBio data of the evolved strains."""
    #Storing all generated dfs per strain for potential future applications.
    genome_length = {strain:None for strain in s.strains}

    #Iterating over all strains
    for strain in s.strains:
        #Get all treatments of a strain
        treatments = s.treatments[strain]
        #df for assembly length
        length = pd.DataFrame(columns=treatments,index=[sample['name'] \
            for sample in s.strains[strain] if sample['platform']== 'pacbio'])
        #df for n contigs
        n_contigs = pd.DataFrame(columns=treatments,index=[sample['name'] \
            for sample in s.strains[strain] if sample['platform']== 'pacbio'])
        
        #Iterating over all samples of a strain
        for sample in s.strains[strain]:
            #Getting genome lenght and n contigs
            if sample['platform'] == 'pacbio':
                #Parsing contigs of assembly
                contigs = [contig for contig in SeqIO.parse(join(sample['dir_name'],'assembly.fasta'),'fasta')]
                #Storing assembly length
                length.at[sample['name'],sample['treatment']] = sum([len(contig) for contig in contigs])
                #Storing n contigs
                n_contigs.at[sample['name'],sample['treatment']] = len(contigs)

        #Wrtiting to dictionary for potential future applications
        genome_length[strain] = length

        #Plotting genome length
        fig = p.subplot_treatments(strain,length)
        #Adding genome length of wild-type
        reference_length = sum([len(contig) for contig in SeqIO.parse(s.references[strain],'fasta')])
        fig.add_hline(y=reference_length,annotation_text='reference',line_dash="dash")
        title = 'Assembly length in '+strain
        #Updating labels and dumping to png
        fig.update_layout(
            xaxis_title='samples',
            yaxis_title='assembly length in bp',
            title=title)
        fig.update_traces(showlegend=False)
        fig.write_image(join('..','plots','genome_length',title.replace(' ','_')+'.png'))
        
        #Plotting n contigs
        fig = p.subplot_treatments(strain,n_contigs)
        #Updating labels and dumping to png
        title = 'N contigs in '+strain
        fig.update_layout(
            xaxis_title='samples',
            yaxis_title='n contigs',
            title=title)
        fig.update_traces(showlegend=False)
        fig.write_image(join('..','plots','contigs',title.replace(' ','_')+'.png'))

    return genome_length

def plot_indel():
    """This function plots the sum of deleted and inserted base pairs."""
    #Storing dfs per strain for potential future applications.
    d = {strain:None for strain in s.strains}
    i = {strain:None for strain in s.strains}

    #Iterating over every strain
    for strain in s.strains:
        #Getting all treatments of a strain
        treatments = s.treatments[strain]
        #df for storing sum of deleted bases
        deleted_bases = pd.DataFrame(columns=treatments,index=[sample['name'] \
            for sample in s.strains[strain] if sample['platform']== 'pacbio'])
        #df for storing sum of inserted bases
        inserted_bases = pd.DataFrame(columns=treatments,index=[sample['name'] \
            for sample in s.strains[strain] if sample['platform']== 'pacbio'])

        #Iterating over every sample of a strain
        for sample in s.strains[strain]:
            if sample['platform'] == 'pacbio':
                #Biggest impact have regions with no coverage outputted from
                #https://github.com/nahanoo/deletion_detection
                no_coverage = join(sample['dir_name'],'no_alignment_regions.tsv')
                if exists(no_coverage):
                    #Writing sum ob base paris with no alignment to df
                    deleted_bases.at[sample['name'],sample['treatment']] = sum(pd.read_csv(no_coverage,sep='\t',\
                        usecols=['chromosome','position','length']).drop_duplicates()['length'])

                in_read = join(sample['dir_name'],'in_read_deletions.tsv')
                if exists(in_read):
                    #Adding sum of base pairs with deletions to df
                    deleted_bases.at[sample['name'],sample['treatment']] += sum(pd.read_csv(in_read,sep='\t',\
                        usecols=['chromosome','position','length']).drop_duplicates()['length'])

                #Storing sum of inserted base pairs to df derived from
                #https://github.com/nahanoo/deletion_detection
                insertions = join(sample['dir_name'],'insertions.tsv')
                if exists(insertions):
                    #Writing sum of inserted base pairs to df
                    inserted_bases.at[sample['name'],sample['treatment']] = sum(pd.read_csv(insertions,sep='\t',\
                        usecols=['chromosome','position','length']).drop_duplicates()['length'])
        
        #Storing dfs in dictionary for future processing
        d[strain] = deleted_bases

        #Plotting deleted base pairs
        fig = p.subplot_treatments(strain,deleted_bases)
        title = 'Deleted bases in '+strain
        fig.update_layout(
            xaxis_title='samples',
            yaxis_title='deleted bp',
            title=title)
        fig.update_traces(showlegend=False)
        fig.write_image(join('..','plots','deleted_bases',title.replace(' ','_')+'.png'))
        
        #Storing dfs in dictionary for future processing
        i[strain] = inserted_bases

        #Plotting inserted base pairs
        fig = p.subplot_treatments(strain,inserted_bases)
        title = 'Inserted bases in '+strain
        fig.update_layout(
            xaxis_title='samples',
            yaxis_title='inserted bp',
            title=title)
        fig.update_traces(showlegend=False)
        fig.write_image(join('..','plots','inserted_bases',title.replace(' ','_')+'.png'))

    return d,i

def get_affected_products(df_name):
    """This function returns all products which were affected
    by deletions or insertions per strain. Counts of observed mutated
    products are summed."""
    affected_products = {strain:None for strain in s.strains}
    #Iterating over every strain
    for strain,samples in s.strains.items():
        #Getting all treatments of a strain
        treatments = s.treatments[strain]
        #Setting up df
        sorted_products = pd.DataFrame(columns=treatments)
        for sample in samples:
            if sample['platform'] == 'pacbio':
                f = join(sample['dir_name'],df_name)
                if exists(f):
                    #Getting all mutated products per sample
                    products = set(pd.read_csv(f,sep='\t').dropna()['product'])
                    for product in products:
                        #Storing how many time product was mutated
                        if product in sorted_products[sample['treatment']].dropna().index:
                            sorted_products.at[product,sample['treatment']] += 1
                        else:
                            sorted_products.at[product,sample['treatment']] = 1
        #Sorting index
        sorted_products.index.name = 'product'
        #Storing as dictionary
        affected_products[strain] = sorted_products

    return affected_products

def deleted_products():
    """This function creates a table listing all products which were
    affected by deletions. It sums the counts of how many times a product
    was mutated."""
    #Getting all deleted products
    affected_products = get_affected_products('no_alignment_regions.tsv')
    #Iterating over every strain
    for strain in s.strains:
        #Getting all treatments per strain
        treatments = s.treatments[strain]
        #Creating full column names
        column_names = ['treatment '+str(treatment) for treatment in treatments]
        fname = 'deleted_products_'+strain.replace(' ','_')+'.csv'
        #Dumping counts of deleted products to csv
        affected_products[strain].to_csv(join('..','tables','pacbio_deletions',fname),header=column_names)

def inserted_products():
    """This function creates a table listing all products which were
    affected by insertions. It sums the counts of how many times a product
    was mutated."""
    #Getting all products affected by insertions
    affected_products = get_affected_products('mutant_to_parent.noalignments.tsv')
    #Iterating over every strain
    for strain in s.strains:
        #Getting all treatments per strain
        treatments = s.treatments[strain]
        #Creating full column names
        column_names = ['treatment '+str(treatment) for treatment in treatments]
        fname = 'inserted_products_'+strain.replace(' ','_')+'.csv'
        #Dumping counts of inserted products to csv
        affected_products[strain].to_csv(join('..','tables','pacbio_insertions',fname),header=column_names)