from numpy import tile
from experiment import Experiment
from samples import Samples
from mutations import Mutations
import plotly.graph_objects as go
import plotly.express as px
from os.path import join

e = Experiment()
s = Samples()
m = Mutations()

class Plotting():
    def plot_gene_count_series(self):
        m.get_gene_count_series()
        styles = ['solid', 'dot', 'dash', 'longdash','solid','dashdot', 'longdashdot','solid']
        for (strain,treatment),df in m.gene_count_series.items():
            fig = go.Figure()
            for i,row in df.iterrows():
                fig.add_trace(go.Scatter(
                    x = e.timepoints,
                    y = row.to_list(),
                    name = i
                ))
            for counter,trace in enumerate(fig.data):
                trace.line.dash = styles[counter]
            fig.update_layout(
                title = ' '.join([strain,'in','treatment',str(treatment)]),
                xaxis_title = 'timepoint',
                yaxis_title = 'observed in n microcosms'
            )
            fig.write_image(join('..','plots','genes','_'.join([strain.replace(' ','_'),str(treatment)])+'.png'))