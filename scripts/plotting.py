from samples import Samples
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd
from os.path import join
from os.path import exists

s = Samples()

class Plotting():
    def subplot_treatments(self,fig,strain,df):
        #Get all treatments of a strain
        treatments = s.treatments[strain]



p = Plotting()

def plot_deletions():
    for strain in ['at','ct']:
        strain = s.abbreviations[strain]
        treatments = s.treatments[strain]
        subplot_titles = ['treatment '+str(treatment) for treatment in treatments]
        fig = make_subplots(rows=1,cols=len(treatments),subplot_titles=subplot_titles)
        out_sum = pd.DataFrame(columns=treatments,index=[sample['name'] \
            for sample in s.strain[strain] if sample['platform']== 'pacbio'])
        for sample in s.strain[strain]':
            if sample['platform'] == 'illumina':
                no_coverage = join(sample['dir_name'],'no_coverage_areas.tsv')
        p.subplot_treatments(fig,strain)