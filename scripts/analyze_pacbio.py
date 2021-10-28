from samples import Samples
from samples import Experiment
from plotting import Plotting
from os.path import join
from os.path import exists
import pandas as pd
from Bio import SeqIO

s = Samples()
p = Plotting()
e = Experiment()

def plot_indel():
    d = {strain:None for strain in s.strains}
    i = {strain:None for strain in s.strains}
    for strain in s.strains:
        treatments = s.treatments[strain]
        deleted_bases = pd.DataFrame(columns=treatments,index=[sample['name'] \
            for sample in s.strains[strain] if sample['platform']== 'pacbio'])
        inserted_bases = pd.DataFrame(columns=treatments,index=[sample['name'] \
            for sample in s.strains[strain] if sample['platform']== 'pacbio'])
        for sample in s.strains[strain]:
            if sample['platform'] == 'pacbio':
                #Summing deleted bases from deletion_detection analysis
                no_coverage = join(sample['dir_name'],'no_alignment_regions.tsv')
                if exists(no_coverage):
                    deleted_bases.at[sample['name'],sample['treatment']] = sum(pd.read_csv(no_coverage,sep='\t',\
                        usecols=['chromosome','position','length']).drop_duplicates()['length'])
                in_read = join(sample['dir_name'],'in_read_deletions.tsv')
                if exists(in_read):
                    deleted_bases.at[sample['name'],sample['treatment']] += sum(pd.read_csv(in_read,sep='\t',\
                        usecols=['chromosome','position','length']).drop_duplicates()['length'])
                insertions = join(sample['dir_name'],'insertions.tsv')
                if exists(insertions):
                    inserted_bases.at[sample['name'],sample['treatment']] = sum(pd.read_csv(insertions,sep='\t',\
                        usecols=['chromosome','position','length']).drop_duplicates()['length'])
        
        d[strain] = deleted_bases
        fig = p.subplot_treatments(strain,deleted_bases)
        title = 'delete bases in '+strain
        fig.update_layout(
            xaxis_title='sample',
            yaxis_title='deleted bp',
            title=title)
        fig.update_traces(showlegend=False)
        fig.write_image(join('..','plots','deleted_bases',title.replace(' ','_')+'.png'))
        
        i[strain] = inserted_bases
        fig = p.subplot_treatments(strain,inserted_bases)
        title = 'inserted bases in '+strain
        fig.update_layout(
            xaxis_title='sample',
            yaxis_title='inserted bp',
            title=title)
        fig.update_traces(showlegend=False)
        fig.write_image(join('..','plots','inserted_bases',title.replace(' ','_')+'.png'))
    return d,i

def plot_genome_length():
    """Plotting assembly length of pacbio data and resulting n contigs of assmblies."""
    genome_length = {strain:None for strain in s.strains}
    for strain in s.strains:
        #Get all treatments of a strain
        treatments = s.treatments[strain]
        #Create subplot titles
        length = pd.DataFrame(columns=treatments,index=[sample['name'] \
            for sample in s.strains[strain] if sample['platform']== 'pacbio'])
        n_contigs = pd.DataFrame(columns=treatments,index=[sample['name'] \
            for sample in s.strains[strain] if sample['platform']== 'pacbio'])
        for counter,treatment in enumerate(treatments):
            #Get all sample names of a treatment
            for sample in s.strains[strain]:
                #Getting genome lenght and n contigs
                if sample['platform'] == 'pacbio':
                    contigs = [contig for contig in SeqIO.parse(join(sample['dir_name'],'assembly.fasta'),'fasta')]
                    n_contigs.at[sample['name'],sample['treatment']] = len(contigs)
                    length.at[sample['name'],sample['treatment']] = sum([len(contig) for contig in contigs])
        
        genome_length[strain] = length
        fig = p.subplot_treatments(strain,length)
        reference_length = sum([len(contig) for contig in SeqIO.parse(s.references[strain],'fasta')])
        fig.add_hline(y=reference_length,annotation_text='reference',line_dash="dash")
        title = 'assembly length in '+strain
        fig.update_layout(
            xaxis_title='sample',
            yaxis_title='assembly length in bp',
            title=title)
        fig.update_traces(showlegend=False)
        fig.write_image(join('..','plots','genome_length',title.replace(' ','_')+'.png'))
        
        fig = p.subplot_treatments(strain,n_contigs)
        title = 'n contigs in '+strain
        fig.update_layout(
            xaxis_title='sample',
            yaxis_title='n contigs',
            title=title)
        fig.update_traces(showlegend=False)
        fig.write_image(join('..','plots','contigs',title.replace(' ','_')+'.png'))
    return genome_length

def plot_differences():
    difference = {strain:None for strain in s.strains}
    d,i = plot_indel()
    g = plot_genome_length()
    for strain in s.strains:
        difference[strain] = g[strain] - d[strain] + i[strain]
        fig = p.subplot_treatments(strain,difference[strain])
        title = 'differences in '+strain
        fig.update_layout(
            xaxis_title='sample',
            yaxis_title='n contigs',
            title=title)
        fig.update_traces(showlegend=False)
        reference_length = sum([len(contig) for contig in SeqIO.parse(s.references[strain],'fasta')])
        fig.add_hline(y=reference_length,annotation_text='reference',line_dash="dash")
        fig.write_image(join('..','plots','differences',title.replace(' ','_')+'.png'))
    return difference