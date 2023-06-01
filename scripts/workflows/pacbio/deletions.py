from samples import Samples
from os.path import join, exists
import pandas as pd
import plotly.express as px
from Bio import SeqIO

s = Samples()

strains = {s.abbreviations['at']: '<i>A. tumefaciens</i>',
           s.abbreviations['ct']: '<i>C. testosteroni</i>',
           s.abbreviations['oa']: '<i>O. anthropi</i>',
           s.abbreviations['ms']: '<i>M. saperdae</i>'}


def parse_deletions():
    df = pd.DataFrame(columns=['strain', 'treatment', 'cosm', 'deletions'])
    i = 0
    species = [s.abbreviations[sp] for sp in ['at','ct']]
    for strain in species:
        for sample in s.strains[strain]:
            if sample['platform'] == 'pacbio':
                ds = 0
                f = join(sample['dir_name'], 'deletions.tsv')
                if exists(f):
                    ds += sum(pd.read_csv(f, sep='\t')['length'])
                f = join(sample['dir_name'], 'plasmids.tsv')
                if exists(f):
                    ds += sum(pd.read_csv(f, sep='\t')['length'])
                df.loc[i] = [strains[strain],
                             sample['treatment'], sample['cosm'], ds]
                i += 1
    return df


def parse_genome_length():
    df = pd.DataFrame(columns=['strain', 'treatment', 'cosm', 'deletions'])
    i = 0
    species = [s.abbreviations[sp] for sp in ['at','ct']]
    for strain in species:
        for sample in s.strains[strain]:
            if sample['platform'] == 'pacbio':
                ds = 0
                f = join(sample['dir_name'], 'assembly.fasta')
                if exists(f):
                    l = sum([len(contig)
                            for contig in SeqIO.parse(f, 'fasta')])
                    df.loc[i] = [strains[strain],
                                 sample['treatment'], sample['cosm'], l]
                    i += 1
    return df

def font_size(fig):
    j = 10
    fig.update_layout(font={'size': j,'color':'black'})
    for a in fig['layout']['annotations']:
        a['font']['size'] = j
        a['font']['color'] = 'black'
    fig['layout']['title']['font']['size'] = j
    fig['layout']['title']['font']['color'] = 'black'
    fig['layout']['legend']['title']['font']['size'] = j
    fig['layout']['legend']['title']['font']['color'] = 'black'
    fig.for_each_xaxis(lambda axis: axis.title.update(font=dict(size=j,color='black')))
    fig.for_each_yaxis(lambda axis: axis.title.update(font=dict(size=j,color='black')))
    for d in fig['data']:
        d['marker']['size'] = 3
        d['line']['width'] = 0.5
    return fig

def plot(df):
    n_colors = len(set(df['treatment']))
    colors = px.colors.sample_colorscale(
        "Agsunset", [n/(n_colors - 1) for n in range(n_colors)])
    fig = px.box(df, x='strain', y='deletions', color='treatment', color_discrete_sequence=colors, log_y=True,
                 points='all', category_orders={'treatment': [1, 2, 3, 4]},height=300,width=350)

    offsetgroups = ['1', '1', '2', '3']
    for i, d in enumerate(fig['data']):
        d['offsetgroup'] = offsetgroups[i]
    fig.update_traces(boxmean=True, quartilemethod="exclusive",
                    pointpos=0, jitter=1)
    fig = font_size(fig)
    return fig


df = parse_deletions()
fig = plot(df)
fig.update_yaxes(title='Deleted base pairs')
fig.update_xaxes(title='Species')
fig.write_image('deletions.svg')


df = parse_genome_length()
fig = plot(df)
fig.update_layout()
fig.update_yaxes(title='Assembly length')
fig.update_xaxes(title='Species')
fig.write_image('length.svg')
