from samples import Samples
from os.path import join, exists
import pandas as pd
import plotly.express as px
from Bio import SeqIO

s = Samples()

strains = {s.abbreviations['at']: 'At',
           s.abbreviations['ct']: 'Ct',
           s.abbreviations['oa']: 'Ms',
           s.abbreviations['ms']: 'Oa'}

# Execute parsing on cluster


def parse_deletions():
    df = pd.DataFrame(columns=['strain', 'treatment', 'cosm', 'deletions'])
    i = 0
    species = [s.abbreviations[sp] for sp in ['at', 'ct']]
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
    df.to_csv(join('..', 'variants', 'deletions.csv'), index=False)


def parse_deletions_illumina():
    df = pd.DataFrame(columns=['strain', 'treatment',
                      'timepoint', 'cosm', 'deletions'])
    species = [s.abbreviations[sp] for sp in ['at', 'ct']]
    i = 0
    for strain in species:
        for sample in s.strains[strain]:
            if sample['platform'] == 'illumina':
                f = join(sample['dir_name'], 'depth.tsv')
                f_marc = join(sample['dir_name'], 'depth_marc.tsv')
                if exists(f_marc):
                    print(sample['name'])
                    ds = pd.read_csv(f_marc, sep='\t', names=[
                                     'contig', 'pos', 'depth'])
                    mask = []
                    for c, d in zip(ds['contig'], ds['depth']):
                        if (c[:2] == s.abbreviations[strain]) & (d == 0):
                            mask.append(True)
                        else:
                            mask.append(False)
                    ds = ds[mask]
                    df.loc[i] = [strains[strain], sample['treatment'],
                                 sample['timepoint'], sample['cosm'], len(ds)]
                    i += 1
    f = join('..', 'variants', 'deletions_illumina.csv')
    f_marc = join('..', 'variants', 'deletions_illumina_marc.csv')
    df.to_csv(f_marc)


def parse_genome_length():
    df = pd.DataFrame(columns=['strain', 'treatment', 'cosm', 'deletions'])
    i = 0
    species = [s.abbreviations[sp] for sp in ['at', 'ct']]
    for strain in species:
        for sample in s.strains[strain]:
            if sample['platform'] == 'pacbio':
                ds = 0
                f = join(sample['dir_name'], 'assembly.contigs.racon.fasta')
                if exists(f):
                    l = sum([len(contig)
                            for contig in SeqIO.parse(f, 'fasta')])
                    df.loc[i] = [strains[strain],
                                 sample['treatment'], sample['cosm'], l]
                    i += 1
    df.to_csv(join('..', 'variants', 'assembly_length.csv'), index=False)

# Local plotting


def font_size(fig):
    j = 10
    fig.update_layout(font={'size': j, 'color': 'black'})
    for a in fig['layout']['annotations']:
        a['font']['size'] = j
        a['font']['color'] = 'black'
    fig['layout']['title']['font']['size'] = j
    fig['layout']['title']['font']['color'] = 'black'
    fig['layout']['legend']['title']['font']['size'] = j
    fig['layout']['legend']['title']['font']['color'] = 'black'
    fig.for_each_xaxis(lambda axis: axis.title.update(
        font=dict(size=j, color='black')))
    fig.for_each_yaxis(lambda axis: axis.title.update(
        font=dict(size=j, color='black')))
    for d in fig['data']:
        d['marker']['size'] = 3
        d['line']['width'] = 0.5
    return fig


def plot(f):
    df = pd.read_csv(f)
    n_colors = len(set(df['treatment']))
    colors = px.colors.sample_colorscale(
        "Agsunset", [n/(n_colors - 1) for n in range(n_colors)])
    fig = px.box(df, x='strain', y='deletions', color='treatment', color_discrete_sequence=colors, log_y=True,
                 points='all', category_orders={'treatment': [1, 2, 3, 4]}, height=300, width=350)

    offsetgroups = ['1', '1', '2', '3']
    for i, d in enumerate(fig['data']):
        d['offsetgroup'] = offsetgroups[i]
    fig.update_traces(boxmean=True, quartilemethod="exclusive",
                      pointpos=0, jitter=1)
    fig = font_size(fig)
    return fig


def plot_deletions(f):
    fig = plot(f)
    fig.update_xaxes(title='')
    fig.update_yaxes(title='Deleted basepairs')
    fig.write_image(join('..', 'plots', 'plots', 'deletions.svg'))
    return fig


def plot_assembly_length(f):
    fig = plot(f)
    fig.update_xaxes(title='')
    fig.update_yaxes(title='Assembly length')
    fig.write_image(join('..', 'plots', 'plots', 'assmblies.svg'))
    return fig


def plot_deletions_illumina(f):
    df = pd.read_csv(f)
    df = df.sort_values(by='treatment', ascending=True)
    n_colors = len(set(df['treatment']))
    colors = px.colors.sample_colorscale(
        "Agsunset", [n/(n_colors - 1) for n in range(n_colors)])
    fig = px.box(df, x='strain', y='deletions', color='treatment', facet_col='timepoint', points='all',
                 category_orders={'timepoint': [
                       'T11', 'T22', 'T33', 'T44']}, color_discrete_sequence=colors,
                 log_y=True, height=250, width=450)
    titles = ['T11', 'T22', 'T33', 'T44']
    for i, t in enumerate(fig['layout']['annotations']):
        t['text'] = titles[i]

    fig.for_each_yaxis(lambda y: y.update(title=''))
    fig['layout']['yaxis']['title']['text'] = 'Deleted basepairs'
    fig.update_xaxes(title=None)
    fig['layout']['legend']['title']['text'] = 'Treatment'
    fig.update_layout(title='', boxgroupgap=0.2, boxgap=0.3)
    fig.update_traces(boxmean=True, quartilemethod="exclusive",
                      pointpos=0, jitter=1)
    fig.for_each_yaxis(lambda yaxis: yaxis.update(rangemode="tozero"))

    offsetgroups = ['1', '1', '1', '1',
                    '1', '1', '1', '1',
                    '2', '2', '2', '2',
                    '3', '3', '3', '3']
    for i, d in enumerate(fig['data']):
        d['offsetgroup'] = offsetgroups[i]

    fig = font_size(fig)
    fig.write_image(join('..','plots','plots','deletions_illumina.svg'))
    return fig



def caller():
    fig = plot_deletions(join('..', 'variants', 'deletions.csv'))
    fig.show()
    fig = plot_assembly_length(join('..', 'variants', 'assembly_length.csv'))
    fig.show()
    fig = plot_deletions_illumina(join('..', 'variants', 'deletions_illumina.csv'))



#fig = plot_deletions_illumina(join('..', 'variants', 'deletions_illumina_marc.csv'))