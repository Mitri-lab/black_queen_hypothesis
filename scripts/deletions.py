from samples import Samples
from os.path import join, exists, split
import pandas as pd
import plotly.express as px
from Bio import SeqIO
import vcfpy
import plotly.graph_objects as go
"""This script hosts the code for the structural anlysis of the PacBio data."""

s = Samples()

strains = {s.abbreviations['at']: 'At',
           s.abbreviations['ct']: 'Ct',
           s.abbreviations['ms']: 'Ml',
           s.abbreviations['oa']: 'Oa'}

colors = {'ct': ['#7570B3', '#E6AB02', '#D95F02'],
          'at': ['#1B9E77', '#E6AB02', '#D95F02'],
          'all': ['#1B9E77', '#7570B3', '#E6AB02', '#D95F02'],
          'cosm': ['#4b2991', '#a431a0', '#ea4f88', '#f89178', '#edd9a3']
          }
colors_t = {'1': '#1B9E77',
            '2': '#7570B3',
            '3': '#E6AB02',
            '4': '#D95F02'}

# Execute parsing on cluster

w = 400
h = 250


def parse_pacbio_coverage():
    """This function stores how many base-pairs don't have coverage.
    """
    df = pd.DataFrame(columns=['strain', 'treatment',
                      'cosm', 'name', 'contig', 'start', 'end', 'length'])
    for strain, samples in s.strains.items():
        for sample in samples:
            if sample['platform'] == 'pacbio':
                f = join(sample['dir_name'], 'depth_Q_0.concat.csv')
                if exists(f):
                    dels = pd.read_csv(f)
                    for i, row in dels.iterrows():
                        df.loc[len(df)] = [strain, sample['treatment'], sample['cosm'], sample['name'], row['chromosome'],
                                           int(row['position']), int(row['position']) + int(row['length']), int(row['length'])]
    df.to_csv(join('..', 'variants', 'pacbio_zero_coverage.csv'))


def parse_pacbio_deletions():
    """This functions stores how many base pairs are inserted or deleted that
    were detected with structural variant calling.
    """
    df = pd.DataFrame(columns=['strain', 'treatment', 'cosm', 'deletions'])
    ins = pd.DataFrame(columns=['strain', 'treatment', 'cosm', 'insertions'])
    for strain, samples in s.strains.items():
        for sample in samples:
            if sample['platform'] == 'pacbio':
                f = join(sample['dir_name'], 'sniffles.var.vcf')
                if exists(f):
                    l = 0
                    k = 0
                    reader = vcfpy.Reader.from_path(f)
                    for record in reader:
                        for i, r in enumerate(record.ALT):
                            if record.INFO['SVTYPE'] == 'DEL':
                                l += record.INFO['SVLEN'] * (-1)
                            if record.INFO['SVTYPE'] == 'INS':
                                k += record.INFO['SVLEN']
                df.loc[len(df)] = [strain, sample['treatment'],
                                   sample['cosm'], l]
                ins.loc[len(ins)] = [
                    strain, sample['treatment'], sample['cosm'], k]
    df.to_csv(join('..', 'variants', 'pacbio_inread_deletions.csv'), index=False)
    ins.to_csv(join('..', 'variants', 'pacbio_inread_insertions.csv'), index=False)


def dump_genes():
    """This function parses how many genes are present in an assembly."""
    df = pd.DataFrame(columns=['name', 'genes', 'treatment', 'strain', 'cosm'])
    # refs = [join(split(r),'bakta','assembly.contigs.polypolish.gbff') for r in s.references

    for strain, samples in s.strains.items():
        for sample in samples:
            if (sample['platform'] == 'pacbio'):
                f = join(sample['dir_name'], 'bakta',
                         'assembly.contigs.racon.tsv')
                if exists(f):
                    tsv = pd.read_csv(f, sep='\t', skiprows=5)
                    tsv = tsv[tsv['Type'] == 'cds']
                    # tsv = tsv.dropna(subset=['Gene'])
                    df.loc[len(df)] = [sample['name'], len(tsv),
                                       sample['treatment'], strain, sample['cosm']]
                else:
                    pass

    df.to_csv(join('..', 'variants', 'number_of_genes.csv'), index=False)



def parse_assemblies():
    """Parses assembly lengths"""
    df = pd.DataFrame(columns=['strain', 'treatment',
                      'cosm', 'assembly_length', 'contigs'])
    i = 0
    species = [s.abbreviations[sp] for sp in ['at', 'ct', 'oa', 'ms']]
    for strain in species:
        for sample in s.strains[strain]:
            if sample['platform'] == 'pacbio':
                ds = 0
                f = join(sample['dir_name'], 'assembly.contigs.racon.fasta')
                if exists(f):
                    contigs = [len(contig)
                               for contig in SeqIO.parse(f, 'fasta')]
                    l, n_contigs = sum(contigs), len(contigs)

                    df.loc[i] = [strains[strain],
                                 sample['treatment'], sample['cosm'], l, n_contigs]
                    i += 1
    df.to_csv(join('..', 'variants', 'assembly_lengths.csv'), index=False)

# Local plotting

def font_size(fig, marker_size=3):
    """Style function for figures setting fot size and true black color."""
    for d in fig['data']:
        d['marker']['size'] = marker_size
        d['line']['width'] = 0.5
    # Font size
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
    fig.update_layout(
        margin=dict(l=20, r=20, t=20, b=20),
    )
    return fig


def plot_insertions_pacbio(abb):
    """Plots how many base pairs were inserted based on structural variant calling."""
    df = pd.read_csv(join('..', 'variants', 'pacbio_inread_insertions.csv'))
    mask = (df['strain'] == s.abbreviations[abb])
    df = df[mask]
    df = df.sort_values(by='treatment')
    df = df.replace(to_replace=0, value=0.1)
    fig = px.box(df, x='treatment', y='insertions', color='treatment',
                 points='all', height=h, width=w/3, log_y=True)
    fig.update_traces(boxmean=True, quartilemethod="linear",
                      pointpos=0, jitter=1)
    fig.update_layout(showlegend=False,
                      title='PacBio inread deletions', title_x=0.5)
    fig.update_xaxes(title='Condition', type='category')
    fig.update_yaxes(title='Deleted bp', rangemode='tozero')
    fig = font_size(fig, marker_size=5)
    fig.update_yaxes(title_standoff=0)
    fig.update_xaxes(title_standoff=0)
    for d in fig['data']:
        d['marker']['color'] = colors_t[d['name']]
        d['line']['color'] = colors_t[d['name']]
    fig.write_image(join('..', 'plots', 'plots', abb +
                    '_pacbio_inread_insertions.svg'))


def plot_deletions_pacbio(abb):
    """Plots how many base pairs were deleted based on structural variant calling."""
    df = pd.read_csv(join('..', 'variants', 'pacbio_inread_deletions.csv'))
    mask = (df['strain'] == s.abbreviations[abb])
    df = df[mask]
    df = df.sort_values(by='treatment')
    # df = df.replace(to_replace=0,value=0.1)
    fig = px.box(df, x='treatment', y='deletions', color='treatment',
                 points='all', height=h, width=w/3, log_y=False)
    fig.update_traces(boxmean=True, quartilemethod="linear",
                      pointpos=0, jitter=1)
    fig.update_layout(showlegend=False, title=abb +
                      ' transfer 44', title_x=0.5)
    fig.update_xaxes(title='Condition', type='category')
    fig.update_yaxes(title='Within read deleted base pairs',
                     rangemode='tozero')
    fig = font_size(fig, marker_size=5)
    fig.update_yaxes(title_standoff=0)
    fig.update_xaxes(title_standoff=0)
    for d in fig['data']:
        d['marker']['color'] = colors_t[d['name']]
        d['line']['color'] = colors_t[d['name']]
    fig.write_image(join('..', 'plots', 'plots', abb +
                    '_pacbio_inread_deletions.svg'))


def plot_assembly_length(abb):
    """
    Figure 3 B Panel C
    Supplementary Figure 11 B Panel C
    Plots assembly lenghts.
    """
    genome_length = {'ct': 6130738,
                     'at': 5430046,
                     'oa': 5168237,
                     'ms': 3761082}


    df = pd.read_csv(join('..', 'variants', 'assembly_lengths.csv'))
    df = df.replace(to_replace='Ct', value='ct')
    df = df.replace(to_replace='At', value='at')
    df = df.replace(to_replace='Oa', value='oa')
    df = df.replace(to_replace='Ml', value='ms')
    mask = (df['strain'] == abb)
    df = df[mask]
    df = df.sort_values(by='treatment')
    # df = df.replace(to_replace=0,value=0.1)
    hover_data = ['treatment', 'cosm']
    fig = px.box(df, x='treatment', y='assembly_length', color='treatment',
                 points='all', height=h, width=w/3, log_y=True, hover_data=hover_data)
    fig.add_hline(genome_length[abb], line_dash="dot",
                  line_color='gray', opacity=0.75)
    fig.update_traces(boxmean=True, quartilemethod="linear",
                      pointpos=0, jitter=1)
    fig.update_layout(showlegend=False, title=abb +
                      ' transfer 44', title_x=0.5)
    fig.update_xaxes(title='Condition', type='category')
    fig.update_yaxes(title='Assebmly length [base pairs]', rangemode='tozero')
    fig = font_size(fig, marker_size=5)
    fig.update_yaxes(title_standoff=0)
    fig.update_xaxes(title_standoff=0)
    for d in fig['data']:
        d['marker']['color'] = colors_t[d['name']]
        d['line']['color'] = colors_t[d['name']]
    fig.write_image(join('..', 'plots', 'plots', abb+'_assembly_length.svg'))
    return fig


def plot_contigs(abb):
    """Plots number of contigs per assembly."""
    genome_length = {'ct': 2,
                     'at': 5,
                     'oa': 4,
                     'ms': 1}

    """ancestors = ['ct','at']
    for a in ancestors:
        l = 0 
        f = join(s.work,'ancestors',a,'assembly.contigs.polypolish.fasta')
        contigs = [contig for contig in SeqIO.parse(f,'fasta')]
        for c in contigs:
            l += len(c)
            genome_length[a] = l"""

    df = pd.read_csv(join('..', 'variants', 'assembly_length.csv'))
    df = df.replace(to_replace='Ct', value='ct')
    df = df.replace(to_replace='At', value='at')
    df = df.replace(to_replace='Oa', value='oa')
    df = df.replace(to_replace='Ml', value='ms')

    mask = (df['strain'] == abb)
    df = df[mask]
    df = df.sort_values(by='treatment')
    # df = df.replace(to_replace=0,value=0.1)
    hover_data = ['treatment', 'cosm']
    fig = px.box(df, x='treatment', y='contigs', color='treatment',
                 points='all', height=h, width=w/3, log_y=True, hover_data=hover_data)
    fig.add_hline(genome_length[abb], line_dash="dot",
                  line_color='gray', opacity=0.75)
    fig.update_traces(boxmean=True, quartilemethod="linear",
                      pointpos=0, jitter=1)
    fig.update_layout(showlegend=False, title=abb +
                      ' transfer 44', title_x=0.5)
    fig.update_xaxes(title='Condition', type='category')
    fig.update_yaxes(title='N contigs', rangemode='tozero')
    fig = font_size(fig, marker_size=5)
    fig.update_yaxes(title_standoff=0)
    fig.update_xaxes(title_standoff=0)
    for d in fig['data']:
        d['marker']['color'] = colors_t[d['name']]
        d['line']['color'] = colors_t[d['name']]
    fig.write_image(join('..', 'plots', 'plots', abb+'_contigs.svg'))


def pacbio_zero_coverage(abb):
    """Supplementary Figure 14 C.
    Plots how many base pairs have zero coverage when mapping the pacbio data
    to the reference genome."""
    df = pd.read_csv(join('..', 'variants', 'pacbio_no_cov.csv'))
    mask = (df['strain'] == s.abbreviations[abb])
    df = df[mask]
    df = df.sort_values(by='treatment')
    # df = df.replace(to_replace=0,value=0.1)
    hover_data = ['treatment', 'cosm']
    fig = px.box(df, x='treatment', y='deletions', color='treatment',
                 points='all', height=h, width=w/3, log_y=True, hover_data=hover_data)
    fig.update_traces(boxmean=True, quartilemethod="linear",
                      pointpos=0, jitter=1)
    fig.update_layout(showlegend=False, title=abb +
                      ' transfer 44', title_x=0.5)
    fig.update_xaxes(title='Condition', type='category')
    fig.update_yaxes(title='Base pairs with 0 PacBio coverage',
                     rangemode='tozero')
    fig = font_size(fig, marker_size=5)
    fig.update_yaxes(title_standoff=0)
    fig.update_xaxes(title_standoff=0)
    for d in fig['data']:
        d['marker']['color'] = colors_t[d['name']]
        d['line']['color'] = colors_t[d['name']]
    fig.write_image(join('..', 'plots', 'plots',
                    abb+'_pacbio_zero_coverage.svg'))
    return fig


def number_of_genes(abb):
    """Supplementary Figure 12 B right panel.
    Plots number of genes in annotated assemblies."""
    df = pd.read_csv(join('..', 'variants', 'genes.csv'))
    mask = (df['strain'] == s.abbreviations[abb])
    df = df[mask]
    df = df.sort_values(by='treatment')
    # df = df.replace(to_replace=0,value=0.1)
    hover_data = ['treatment', 'cosm']
    fig = px.box(df, x='treatment', y='genes', color='treatment',
                 points='all', height=h, width=w/3, log_y=True, hover_data=hover_data)
    fig.update_traces(boxmean=True, quartilemethod="linear",
                      pointpos=0, jitter=1)
    fig.update_layout(showlegend=False, title=abb +
                      ' transfer 44', title_x=0.5)
    fig.update_xaxes(title='Condition', type='category')
    fig.update_yaxes(title='Number of coding sequences', rangemode='tozero')
    fig = font_size(fig, marker_size=5)
    fig.update_yaxes(title_standoff=0)
    fig.update_xaxes(title_standoff=0)
    for d in fig['data']:
        d['marker']['color'] = colors_t[d['name']]
        d['line']['color'] = colors_t[d['name']]
    fig.write_image(join('..', 'plots', 'plots', abb+'_genes.svg'))



def large_deletion():
    """Supplementary Figure S14 C
    Plots the coverage for the two samples with the large deletion.
    Only works on the cluster
    """
    dfs = []
    df = pd.read_csv(join(s.work, 'Ct45.1', 'depth_Q_0.tsv'), sep='\t')
    df.columns = ['chromosome', 'position', 'coverage']
    df.insert(len(df.columns), 'sample', 'Ct45.1')
    dfs.append(df)

    df = pd.read_csv(join(s.work, 'Ct42.2', 'depth_Q_0.tsv'), sep='\t')
    df.columns = ['chromosome', 'position', 'coverage']
    df.insert(len(df.columns), 'sample', 'Ct42.2')
    dfs.append(df)

    df = pd.concat(dfs)
    df = df[df['chromosome'] == 'tig00000001_polypolish']

    fig = px.scatter(df, x='position', y='coverage',
                     color='sample', width=w, height=h, opacity=0.5)
    fig = font_size(fig)
    fig.write_image(join('..','plots','plots','pacbio_deletion.svg'))
