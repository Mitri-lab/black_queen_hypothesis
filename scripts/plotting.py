from os.path import join
from samples import Samples
from samples import Experiment
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px
import random

"""
################################################################################
Author: https://github.com/nahanoo
This script acts as a small plotting library used to visualize findings in
PacBio and Illumina data. Most of the time we look at the data from a treatment
perspective, which is why most of the plots are generated using
Plotting().subplot_treatments. The plots are generated in analyze_illumina.py
and analyze_pacbio.py. The plots can be inspected in ../plots are or in
../reports/report.pdf if you also want some context.
################################################################################
"""

s = Samples()
color = 'mediumpurple'
labels = {'at': '<i>A. tumefaciens</i>',
              'ct': '<i>C. testosteroni</i>',
              'ms': '<i>M. saperdae</i>',
              'oa': '<i>O. anthropi</i>'}


class Plotting():
    """This class contains all the plotting backbone for this analysis."""

    def subplot_strains(self,d):
        fig = go.Figure()
        for counter, strain in enumerate(s.strains):
            df = d[strain]
            xs = list(df.columns) * len(df)
            ys = []
            for row in df.values:
                for item in row:
                    ys.append(item)
            name = labels[s.abbreviations[strain]]
            fig.add_trace(go.Scatter(x=xs, y=ys,
                        mode='markers',name=name))
        fig.update_xaxes(type='category')
        fig.update_xaxes(categoryorder='array', categoryarray= [1,2,3,4])
        return fig
        
    def subplot_violin(self,d):
        titles = [labels[s.abbreviations[strain]] for strain in s.strains]
        fig = make_subplots(rows=1, cols=len(s.strains), shared_yaxes=False,
                                subplot_titles=titles)
        
        for counter, strain in enumerate(s.strains):
            df = d[strain]
            xs = list(df.columns) * len(df)
            ys = []
            for row in df.values:
                for item in row:
                    ys.append(item)

            fig.add_trace(go.Violin(x=xs, y=ys,
                        points='all',pointpos=0,marker_color=color,spanmode='hard'), row=1, col=counter+1)
        fig.update_layout(
            xaxis_title='Treatments',
            yaxis_title='SNPs ',
            margin=dict(
                l=0,
                r=10,
                b=0,
                t=45,
                pad=4
                ),
            width=800
            )
        #fig.update_yaxes(type='log')
        fig.update_xaxes(type='category')
        fig.update_traces(showlegend=False)
        fig.update_traces({'fillcolor': 'rgba(255,255,255,0)',
        'line': {'color': 'rgba(255,255,255,0)'}})
        fig.write_image(join('..', 'plots', 'snps', 'snps_all2.png'), scale=2)

        return fig

    def plot_effects(self,snps):
        label_added = []

        colors = {'synonymous_variant':'#1f77b4',
            'missense_variant':'#7f7f7f',
            'disruptive_inframe_deletion':'#9467bd',
            'stop_gained':'#e377c2',
            'conservative_inframe_deletion':'#17becf',
            'conservative_inframe_insertion':'blue',
            'frameshift_variant':'purple'}

        for counter, strain in enumerate(s.strains):
            fig = go.Figure()
            df = snps[strain]
            for effect,row in df.iterrows():
                xs = list(df.columns)
                ys = row.to_list()
                try:
                    fig.add_trace(go.Violin(x=xs,y=ys,name=effect.replace('_',' '),\
                                pointpos=0,points='all',spanmode='hard', legendgroup=effect,marker_color=colors[effect]))
                except KeyError:
                    pass
            fig.update_layout(
                xaxis_title='Treatments',
                yaxis_title='SNPs ',
                margin=dict(
                        l=0,
                        r=10,
                        b=0,
                        t=45,
                        pad=4
                        ),
                    width=400,
                    height=270
                )
            #fig.update_yaxes(type='log')
            fig.update_layout(violinmode='group')
            fig.update_xaxes(type='category')
            fig.update_traces({'fillcolor': 'rgba(255,255,255,0)',
                    'line': {'color': 'rgba(255,255,255,0)'}})
            fig.write_image(join('..', 'plots', 'snps',s.abbreviations[strain] + 'snps.png'), scale=2)

    def subplot_treatments(self,strain,df):
        """This function creates subplots where one subplot is one treatment.
        On the x-axis all samples are listed with the name stored in sample['name'].
        Input is a df with treatments as columns and sample names as index."""
        #Get all treatments of a strain
        treatments = s.treatments[strain]
        #Create subplot titles
        subplot_titles = ['treatment '+str(treatment) for treatment in treatments]
        #Subplot with shared yaxes
        fig = make_subplots(rows=1,cols=len(treatments),subplot_titles=subplot_titles,shared_yaxes=True)
        #Adding every treatment as trace
        for counter,treatment in enumerate(treatments):
            #Dropping nas
            subset = df[treatment].dropna()
            fig.add_trace(go.Scatter(x=subset.index,y=subset.values,mode='markers',marker_color=color),\
                row=1,col=counter+1)
        return fig

    def subplot_products(self,strain,df):
        """This function plots the products impacted by mutations.
        It creates subplots where one suplot is one treatment.
        Input is a df with treatments as columns and products as index.
        Values are observed counts."""
        #Get all treatments of a strain
        treatments = s.treatments[strain]
        subplot_titles = ['treatment '+str(treatment) for treatment in treatments]
        #N products varies a lot which is why we adjust row width.
        height_fractions = [len(df[treatment].dropna())/len(df) for treatment in treatments]
        fig = make_subplots(rows=len(treatments),cols=1,subplot_titles=subplot_titles,\
            shared_xaxes=False,row_heights=height_fractions)
        #Adding every treatment as trace
        for counter,treatment in enumerate(treatments):
            #Dropping na and sorting index
            subset = df[treatment].dropna().sort_index()
            fig.add_trace(go.Scatter(x=subset.values,y=subset.index,mode='markers',marker_color=color),\
                row=counter+1,col=1)
        fig.update_yaxes(automargin=True)
        return fig

    def trajectories(self,df):
        """This function plots mutated gene counts over
        every timepoint and microcosms. It's presented as a timeline.
        Input is a df with timepoints as columns and gene names as index"""
        #Grabbing timepoints from experiment class
        e = Experiment()
        #Sometimes we have many genes which is why we need different styles
        styles = ['dot', 'dash', 'longdash','dashdot', 'longdashdot','solid']
        fig = go.Figure()
        df = df.sort_index()
        #Adding gene counts for every gene as trace
        for i,row in df.iterrows():
            fig.add_trace(go.Scatter(
                x = e.timepoints,
                y = row.to_list(),
                name = i)
            )
        #Picking random style for every gene
        for counter,trace in enumerate(fig.data):
            trace.line.dash = random.choice(styles)
        return fig

    