from samples import Samples
from samples import Experiment
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import random

s = Samples()
color = 'mediumpurple'

class Plotting():
    """This class contains all the plotting structure for this analysis."""   
    def subplot_treatments(self,strain,df):
        """This function creates a subplot with a subplot per treatment.
        Input is a df with treatments as columns and sample names as index."""
        #Get all treatments of a strain
        treatments = s.treatments[strain]
        subplot_titles = ['treatment '+str(treatment) for treatment in treatments]
        #We create the figure here because we use subplot_treatments also with shared axes
        fig = make_subplots(rows=1,cols=len(treatments),subplot_titles=subplot_titles,shared_yaxes=True)
        for counter,treatment in enumerate(treatments):
            subset = df[treatment].dropna()
            fig.add_trace(go.Scatter(x=subset.index,y=subset.values,mode='markers',marker_color=color),\
                row=1,col=counter+1)
        return fig

    def subplot_products(self,strain,df):
        """This function plots the products impacted by mutations.
        Input is a df with treatments as columns and products as index.
        Values are observed counts."""
        #Get all treatments of a strain
        treatments = s.treatments[strain]
        subplot_titles = ['treatment '+str(treatment) for treatment in treatments]
        fig = make_subplots(rows=len(treatments),cols=1,subplot_titles=subplot_titles,\
            shared_xaxes=False,vertical_spacing=0.1)
        for counter,treatment in enumerate(treatments):
            subset = df[treatment].dropna()
            fig.add_trace(go.Scatter(x=subset.values,y=subset.index,mode='markers',marker_color=color),\
                row=counter+1,col=1)
        fig.update_yaxes(automargin=True)
        return fig

    def trajectories(self,df):
        """This function plots tracetories of mutated genes across microcosms.
        Input is a df with timepoints as columns and gene names as index"""
        styles = ['dot', 'dash', 'longdash','dashdot', 'longdashdot','solid']
        fig = go.Figure()
        e = Experiment()
        for i,row in df.iterrows():
            fig.add_trace(go.Scatter(
                x = e.timepoints,
                y = row.to_list(),
                name = i)
            )
            for counter,trace in enumerate(fig.data):
                trace.line.dash = random.choice(styles)
        return fig