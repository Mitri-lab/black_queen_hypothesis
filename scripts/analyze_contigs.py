from Bio import SeqIO
from samples import Samples
from os.path import join
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import pandas as pd

s = Samples()
work = '/work/FAC/FBM/DMF/smitri/evomicrocomm/genome_size/data/'

def plot_contigs(strain):
    """Plots genome length and n contigs per treatment"""
    #Get all treatments of a strain
    treatments = s.treatments[strain]
    #Create subplot titles
    subplot_titles = ['treatment '+str(treatment) for treatment in treatments]
    #Creating a subplot fer every treatment
    fig = make_subplots(rows=1,cols=3,subplot_titles=subplot_titles,\
        specs=[[{"secondary_y": True},{"secondary_y": True},{"secondary_y": True}]],shared_yaxes=True,)
    #Hacky approach to show legend only once
    legend_shown = False
    for counter,treatment in enumerate(treatments):
        #Get all sample names of a treatment
        names = [sample['name'] for sample in s.strains[strain] \
            if (sample['platform'] == 'pacbio') & (sample['treatment'] == treatment)]
        df = pd.DataFrame(columns=['length','contigs'],index=names)
        for sample in s.strains[strain]:
            #Getting genome lenght and n contigs
            if (sample['platform'] == 'pacbio') & (sample['treatment'] == treatment):
                l = 0
                contigs = [contig for contig in SeqIO.parse(join(sample['dir_name'],'assembly.fasta'),'fasta')]
                df.at[sample['name'],'contigs'] = len(contigs)
                for contig in contigs:
                    l += len(contig)
                df.at[sample['name'],'length'] = l
        #Showling legend only once
        if not legend_shown:
            #Adding scatter trace for genome length
            fig.add_trace(go.Scatter(x=df.index,y=df['length'],mode='markers',\
                name='length',opacity=0.5,showlegend=True,marker_color='blue'),row=1,col=counter+1,secondary_y=False)
            #Adding bar trace for n contigs
            fig.add_trace(go.Bar(x=df.index,y=df['contigs'],name='contigs',opacity=0.5,showlegend=True,\
                marker_color='mediumpurple'),row=1,col=counter+1,secondary_y=True)
            legend_shown=True
        if legend_shown:
            fig.add_trace(go.Scatter(x=df.index,y=df['length'],mode='markers',\
                name='length',opacity=0.5,showlegend=False,marker_color='blue'),row=1,col=counter+1,secondary_y=False)
            fig.add_trace(go.Bar(x=df.index,y=df['contigs'],name='contigs',opacity=0.5,showlegend=False,\
                marker_color='mediumpurple'),row=1,col=counter+1,secondary_y=True)

        reference_length = sum([len(contig) for contig in SeqIO.parse(join(work,s.abbreviations[strain],'reference.fasta'),'fasta')])

        fig.add_hline(y=reference_length,annotation_text='reference',line_dash="dash")
        #Adding axixs titles
        fig.update_yaxes(title_text='genome length in bp',secondary_y=False)
        fig.update_yaxes(title_text='n contigs',secondary_y=True)
        yaxis = ['yaxis'+str(n) for n in range(2,6)]
        for axis in yaxis:
            fig['layout'][axis]['title']['text'] = ''
    fig.write_image(join('..','plots','contigs',strain+'.png'))
    
for strain in ['at','ct']:
    fig = plot_contigs(s.abbreviations[strain])