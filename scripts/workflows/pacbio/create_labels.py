import json
import sys
from os.path import join

def create_plot_labels(sample,strain,out):
    """Creating labels and storing as json for
    https://github.com/nahanoo/gc_bias/tree/main/gc_bias
    """
    labels = dict()
    density_plot = dict()
    density_plot['title'] = strain+' in '+sample
    density_plot['xlabel'] = 'coverage per window'
    density_plot['ylabel'] = 'GC content per window'
    density_plot['theme'] = 'plotly'
    labels['density_plot'] = density_plot

    histogram = dict()
    histogram['title'] = strain+' in '+sample
    histogram['xlabel'] = 'coverage'
    histogram['ylabel'] = 'counts'
    histogram['theme'] = 'plotly'
    labels['histogram'] = histogram

    j = json.dumps(labels,indent=4) 
    with open(join(out,'labels.json'),'w') as handle:
        handle.write(j)

if __name__ == "__main__":
    create_plot_labels(sys.argv[1],' '.join([sys.argv[2],sys.argv[3]]),sys.argv[4])
