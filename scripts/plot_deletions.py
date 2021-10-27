import pandas as pd
from samples import Samples
import os
import plotly.express as px
import plotly.graph_objects as go
import json
import subprocess

class Plotting():
    """Plotting class for black queen hypothesis analysis."""
    def categorical_plot(self,df,):
        """This dot plot function takes a df as argument which typically has
        treatments as columns and samples as rows."""
        x = []
        y = []
        for i,row in df.iterrows():
            for counter,value in enumerate(row):
                if not pd.isna(value):
                    x.append(df.columns[counter])
                    y.append(value)
        #First number of size arguments is 16:9 ratio
        self.strip_plot = px.strip(x=x,y=y,width=1057.92*0.8*0.6,height=595.2*0.8)
        self.strip_plot.update_xaxes(type='category')

    def update_labels(self,plot,title,xlabel,ylabel,legend_title):
        """This function can be used to update all labels of a plotly figure."""
        plot.update_layout(
            title=title,
            xaxis_title=xlabel,
            yaxis_title=ylabel,
            legend_title=legend_title)
        plot.update_traces(showlegend=False)


def parse_deletions(meta_dict):
    """Parses all deletions from PacBio analysis."""
    s = Samples()
    base_dir = '/work/FAC/FBM/DMF/smitri/evomicrocomm/genome_size/data/'
    #Here we get all the treatments of one strain
    treatments = list(set([sample['treatment'] for sample in meta_dict \
        if sample['platform']== 'pacbio']))
    #Creating out dataframe for deleted basepairs
    out_sum = pd.DataFrame(columns=treatments,index=[sample['name'] for sample in meta_dict \
        if sample['platform']== 'pacbio'])
    #Creating out datafram for n deletions
    out_n = pd.DataFrame(columns=treatments,index=[sample['name'] for sample in meta_dict \
        if sample['platform']== 'pacbio'])
    for sample in meta_dict:
        #Analysis is only for pacbio data
        if sample['platform'] == 'pacbio':
            #Getting dfs for no aligments as well as inread deletions
            in_read_deletions = pd.read_csv(os.path.join(base_dir,sample['dir_name'],\
                'in_read_deletions.tsv'),sep='\t',usecols=['chromosome','position','length']).drop_duplicates()
            no_alignments = pd.read_csv(os.path.join(base_dir,sample['dir_name'],\
                'no_alignment_regions.tsv'),sep='\t',usecols=['chromosome','position','length']).drop_duplicates()
            df = pd.concat([in_read_deletions,no_alignments])
            out_sum.at[sample['name'],sample['treatment']] = df['length'].sum()
            out_n.at[sample['name'],sample['treatment']] = len(df)
    return out_sum,out_n

def plot_deletions():
    """Plotting pacbio deletions."""
    p = Plotting()
    s = Samples()

    for strain in s.strains:
        #Getting dfs per strain regarding deleted bases and n deletions
        out_sum,out_n = parse_deletions(s.strains[strain])
        #Plotting deleted bases
        p.categorical_plot(out_sum)
        p.update_labels(p.strip_plot,'Deleted bases '+strain,'treatment','deleted basepairs','evolved samples')
        fname = (strain+'_deleted bases'+'.png').replace(' ','_')
        p.strip_plot.write_image(os.path.join('../plots','deleted_bases',fname))
        
        #Plotting n deletions
        p.categorical_plot(out_n)
        p.update_labels(p.strip_plot,'Number of deletions '+strain,'treatment','number of deletions','evolved samples')
        fname = (strain+'_number of deletions'+'.png').replace(' ','_')
        p.strip_plot.write_image(os.path.join('../plots','deleted_bases',fname))
