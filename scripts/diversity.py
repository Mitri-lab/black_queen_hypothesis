from scipy.stats import ttest_ind, kruskal
from samples import Samples
from os.path import join, exists
import vcfpy
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import math
from scipy.spatial import distance
from plotly.subplots import make_subplots
# Sample class for parsing HPC directories contact me for infos
# eric.ulrich@unil.ch
s = Samples()

# Global names for species used for plotting
strains = {s.abbreviations['at']: 'At',
           s.abbreviations['ct']: 'Ct',
           s.abbreviations['oa']: 'Oa',
           s.abbreviations['ms']: 'Ms'}

# Species colors
colors = {'ct': ['#7570B3', '#E6AB02', '#D95F02'],
          'at': ['#1B9E77', '#E6AB02', '#D95F02'],
          'ms': ['#E6AB02', '#D95F02'],
          'oa': ['#D95F02'],
          'all': ['#1B9E77', '#7570B3', '#E6AB02', '#D95F02'],
          'cosm': ['#4b2991', '#a431a0', '#ea4f88', '#f89178', '#edd9a3']
          }

# Condition colors
colors_t = {'1': '#1B9E77',
            '2': '#7570B3',
            '3': '#E6AB02',
            '4': '#D95F02'}

# Plot sizes
w = 400
h = 250


def font_size(fig, marker_size=3, line_size=0.5, fsize=10):
    """Style function for figures setting fot size and true black color."""
    for d in fig['data']:
        d['marker']['size'] = marker_size
        d['line']['width'] = line_size
    # Font size
    j = fsize
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

# This section is executed on the cluster for parsing data on cluster
# To check how the vcf files were generated, check the Snakemake file /workflows/illumina/Snakemake


def parse_variants():
    """Parses vcf files from Illumina data and dumps variants with metdadata as csv."""
    dfs = []
    for strain, samples in s.strains.items():
        for sample in samples:
            if sample['platform'] == 'illumina':
                f = join(sample['dir_name'], 'var.annot.vcf')
                if exists(f):
                    reader = vcfpy.Reader.from_path(f)
                    for record in reader:
                        # Iterating over each variant in case different variants
                        # are present at the same position
                        for i, r in enumerate(record.ALT):
                            snp = {}
                            # Chromosome
                            snp['chrom'] = record.CHROM
                            # Position
                            snp['pos'] = record.POS
                            # Variant quality
                            snp['qual'] = record.QUAL
                            # Sequencing depth
                            snp['depth'] = record.INFO['DP']
                            # Variant frequency
                            snp['freq'] = record.INFO['AO'][i] / \
                                record.INFO['DP']
                            # Combined frequency
                            snp['total_freq'] = sum(
                                record.INFO['AO']) / record.INFO['DP']
                            # Alt sequence
                            snp['alt'] = r
                            # Alt depth
                            snp['alt_count'] = record.INFO['AO'][i]
                            # Reference sequence
                            snp['ref'] = record.REF
                            # Mutation type
                            snp['type'] = record.INFO['TYPE'][i]
                            # Mutation length for INS and DEL
                            snp['len'] = record.INFO['LEN'][i]
                            # Coding or non-coding
                            snp['eff'] = record.INFO['ANN'][i].split('|')[1]
                            # Creates unique key per variant used for plotting trajectories
                            key = '.'.join([snp['chrom'], str(snp['pos']),
                                            str(sample['treatment']), str(sample['cosm']), str(snp['alt'])])
                            snp['linegroup'] = key
                            # Quality cutoff
                            if snp['qual'] >= 20:
                                df = (pd.concat(
                                    [pd.DataFrame(snp, index=[0]), pd.DataFrame(sample, index=[0])], axis=1))
                                dfs.append(df)

    out = pd.concat(dfs)
    out['treatment'] = out['treatment'].astype(str)
    out['cosm'] = out['cosm'].astype(str)
    # Here we apply different filters
    # Variants
    # Non synonymous variants plus alt_count and freq cut offs
    filter = (out['alt_count'] >= 3) & (out['freq'] >= 0.1) & (
        out['qual'] >= 20) & (out['eff'] != 'synonymous_variant')
    out[filter].to_csv(
        join('..', 'variants', 'ns_variants_illumina.csv'), index=False)

    # Synonymous and non synonymous variants
    filter = (out['alt_count'] >= 3) & (out['freq'] >= 0.1) & (
        out['qual'] >= 20)
    out[filter].to_csv(
        join('..', 'variants', 'variants_illumina.csv'), index=False)

    # Fixed variants non synonymous
    filter = (out['alt_count'] >= 3) & (out['total_freq'] >= 0.95) & (
        out['qual'] >= 20) & (out['eff'] != 'synonymous_variant')
    out[filter].to_csv(
        join('..', 'variants', 'fixed_ns_variants_illumina.csv'), index=False)

    # Fixed variants synonymous and non synonymous
    filter = (out['alt_count'] >= 3) & (out['total_freq'] >= 0.95) & (
        out['qual'] >= 20)
    out[filter].to_csv(
        join('..', 'variants', 'fixed_variants_illumina.csv'), index=False)


def parse_snps():
    """Parses fixed SNPs from snippy. Data from snippy is not used for the paper.
    For comments check parse_variants"""
    dfs = []
    for strain, samples in s.strains.items():
        for sample in samples:
            if sample['platform'] == 'illumina':
                f = join(sample['dir_name'], 'snippy', 'snps.vcf')
                if exists(f):
                    reader = vcfpy.Reader.from_path(f)
                    for record in reader:
                        for i, r in enumerate(record.ALT):
                            snp = {}
                            snp['chrom'] = record.CHROM
                            snp['pos'] = record.POS
                            snp['qual'] = record.QUAL
                            snp['depth'] = record.INFO['DP']
                            snp['freq'] = record.INFO['AO'][i] / \
                                record.INFO['DP']
                            snp['alt'] = r
                            snp['ref'] = record.REF
                            key = '.'.join([snp['chrom'], str(snp['pos']),
                                            str(sample['treatment']), str(sample['cosm']), str(snp['alt'])])
                            snp['linegroup'] = key
                            snp['type'] = record.INFO['TYPE'][i]
                            snp['eff'] = record.INFO['ANN'][i].split('|')[1]
                            if snp['qual'] >= 20:
                                df = (pd.concat(
                                    [pd.DataFrame(snp, index=[0]), pd.DataFrame(sample, index=[0])], axis=1))
                                dfs.append(df)
    out = pd.concat(dfs)
    out = out[out['freq'] != 0]
    out['treatment'] = out['treatment'].astype(str)
    out['cosm'] = out['cosm'].astype(str)
    out.to_csv(join('..', 'variants', 'snps_illumina'), index=False)
    out = out[out['eff'] != 'synonymous_variant']
    out.to_csv(join('..', 'variants', 'ns_snps_illumina.csv'), index=False)


def parse_snps_pacbio():
    """Parses SNPs from pacbio data.
    Clonal data, SNPs identified with assemblies.
    """
    dfs = []
    for strain, samples in s.strains.items():
        for sample in samples:
            if sample['platform'] == 'pacbio':
                f = join(sample['dir_name'], 'snippy', 'snps.vcf')
                if exists(f):
                    reader = vcfpy.Reader.from_path(f)
                    for record in reader:
                        for i, r in enumerate(record.ALT):
                            snp = {}
                            snp['chrom'] = record.CHROM
                            snp['pos'] = record.POS
                            snp['qual'] = record.QUAL
                            snp['depth'] = record.INFO['DP']
                            snp['freq'] = record.INFO['AO'][i] / \
                                record.INFO['DP']
                            snp['alt'] = r
                            key = '.'.join([snp['chrom'], str(snp['pos']),
                                            str(sample['treatment']), str(sample['cosm']), str(snp['alt'])])
                            snp['linegroup'] = key
                            if snp['qual'] >= 20:
                                if (sample['strain'] == s.abbreviations['oa']) & (snp['chrom'] == 'tig00000002_polypolish') & (snp['pos'] in range(119960, 120700)):
                                    pass
                                else:
                                    df = (pd.concat(
                                        [pd.DataFrame(snp, index=[0]), pd.DataFrame(sample, index=[0])], axis=1))
                                    dfs.append(df)
    out = pd.concat(dfs)
    out = out[out['freq'] != 0]
    out['treatment'] = out['treatment'].astype(str)
    out['cosm'] = out['cosm'].astype(str)
    out.to_csv(join('..', 'variants', 'snps_pacbio.csv'), index=False)


def depth():
    """Calculates mean coverage using the output 
    from samtools depth -aa Q0 (check Snakemake for details)."""
    out = pd.DataFrame(columns=['sample', 'timepoint',
                                'species', 'treatment', 'depth', 'cosm'])
    i = 0
    for strain, samples in s.strains.items():
        for sample in samples:
            if sample['platform'] == 'illumina':
                f = join(sample['dir_name'], 'depth_Q_0.tsv')
                if exists(f):
                    df = pd.read_csv(f, sep='\t')
                    df.columns = ['chromosome', 'position', 'depth']
                    mean_cov = sum(df['depth']) / len(df)
                    row = [sample['name'], sample['timepoint'], strains[strain],
                           sample['treatment'], mean_cov, sample['cosm']]
                    out.loc[len(out)] = row

    out.to_csv(join('..', 'variants', 'mean_coverage.csv'), index=False)


def mean_coverage():
    """
    Supplementary figure 12
    Plots mean illumina coverage for every species and every timepoint
    """
    df = pd.read_csv(join('..', 'variants', 'mean_coverage.csv'))
    df = df.sort_values(by='treatment', ascending=True)
    # In case of depth 0 remove from data frame for log scale
    depth = [1 if i <= 1 else i for i in df['depth']]
    df['depth'] = depth
    hover_data = ['treatment', 'timepoint', 'species', 'depth', 'cosm']
    fig = px.box(df, x='species', y='depth', color='treatment', facet_col='timepoint', points='all',
                 category_orders={'timepoint': [
                       'T11', 'T22', 'T33', 'T44']}, color_discrete_sequence=colors['all'],
                 log_y=True, hover_data=hover_data, height=h, width=w*2)
    titles = ['Transfer 11', 'Transfer 22', 'Transfer 33', 'Transfer 44']
    for i, t in enumerate(fig['layout']['annotations']):
        t['text'] = titles[i]
    fig.for_each_yaxis(lambda y: y.update(title=''))
    fig['layout']['yaxis']['title']['text'] = 'Mean coverage'
    fig.update_xaxes(title=None)
    fig['layout']['legend']['title']['text'] = 'Condition'
    fig.update_layout(title='', boxgroupgap=0.2, boxgap=0.3)
    fig.update_traces(boxmean=True, quartilemethod="linear",
                      pointpos=0, jitter=1)
    # fig.for_each_yaxis(lambda yaxis: yaxis.update(rangemode="tozero"))

    offsetgroups = ['1', '1', '1', '1',
                    '1', '1', '1', '1',
                    '2', '2', '2', '2',
                    '3', '3', '3', '3']
    for i, d in enumerate(fig['data']):
        d['offsetgroup'] = offsetgroups[i]

    fig = font_size(fig)
    fig.update_yaxes(title_standoff=0)
    fig.update_xaxes(title_standoff=0)
    fig.for_each_yaxis(lambda yaxis: yaxis.update(rangemode="tozero"))
    fig.write_image(join('..', 'plots', 'plots',
                    'mean_coverage.svg'))
    return fig


def n_variants(abb):
    """
    Figure 3 B Panel A
    Figure S11 B Panel A
    Plots the number of variants per microcosm for a species at the last timepoint.
    """
    df = pd.read_csv(join('..', 'variants', 'variants_illumina.csv'))
    mask = (df['strain'] == s.abbreviations[abb]) & (df['timepoint'] == 'T44')
    df = df[mask]
    out = pd.DataFrame(columns=['strain', 'treatment', 'n_variants', 'cosm'])
    for sample in s.strains[s.abbreviations[abb]]:
        if (sample['platform'] == 'illumina') & (sample['timepoint'] == 'T44'):
            n_variants = len(df[df['name'] == sample['name']])
            out.loc[len(out)] = [s.abbreviations[abb],
                                 sample['treatment'], n_variants, sample['cosm']]
    fig = px.box(out, x='treatment', y='n_variants', color='treatment', color_discrete_sequence=colors[abb],
                 points='all', category_orders={'timepoint': ['T11', 'T22', 'T33', 'T44']}, height=h, width=w/3)
    fig.update_traces(boxmean=True, quartilemethod="linear",
                      pointpos=0, jitter=1)
    fig.update_layout(
        showlegend=False, title=strains[s.abbreviations[abb]] + ' transfer 44', title_x=0.5)
    fig.update_xaxes(title='Condition', type='category')
    fig.update_yaxes(title='Number of variants', rangemode='tozero')
    fig = font_size(fig, marker_size=5)
    fig.update_yaxes(title_standoff=0)
    fig.update_xaxes(title_standoff=0)
    fig.write_image(join('..', 'plots', 'plots',
                    abb+'_number_of_variants.svg'))

def ratio_fixed_variants(abb):
    """
    Figure 3 B Panel B
    Figure S11 B Panel B
    Plots the proportion of fixed variants.
    """
    df = pd.read_csv(join('..', 'variants', 'variants_illumina.csv'))
    abb = 'ct'
    mask = (df['strain'] == s.abbreviations[abb]) & (df['timepoint'] == 'T44')
    variants = df[mask]
    df = pd.read_csv(join('..', 'variants', 'fixed_variants_illumina.csv'))
    abb = 'ct'
    mask = (df['strain'] == s.abbreviations[abb]) & (df['timepoint'] == 'T44')
    snps = df[mask]
    out = pd.DataFrame(columns=['strain', 'treatment', 'proportion', 'cosm'])
    for sample in s.strains[s.abbreviations[abb]]:
        if (sample['platform'] == 'illumina') & (sample['timepoint'] == 'T44'):
            proportion = len(snps[snps['name'] == sample['name']]) / \
                len(variants[variants['name'] == sample['name']])
            out.loc[len(out)] = [s.abbreviations[abb],
                                sample['treatment'], proportion, sample['cosm']]
    fig = px.box(out, x='treatment', y='proportion', color='treatment', color_discrete_sequence=colors[abb],
                    points='all', category_orders={'timepoint': ['T11', 'T22', 'T33', 'T44']}, height=h, width=w/3)
    fig.update_traces(boxmean=True, quartilemethod="linear",
                        pointpos=0, jitter=1)
    fig.update_layout(
        showlegend=False, title=strains[s.abbreviations[abb]] + ' transfer 44', title_x=0.5)
    fig.update_xaxes(title='Condition', type='category')
    fig.update_yaxes(title='Proportion of fixed variants', rangemode='tozero')
    fig = font_size(fig, marker_size=5)
    fig.update_yaxes(title_standoff=0)
    fig.update_xaxes(title_standoff=0)
    fig.write_image(join('..', 'plots', 'plots', abb+'_proportion.svg'))

def box_treatments(abb):
    df = get_mutations(add_T0=False)
    filter = (df['strain'] == strains[s.abbreviations[abb]]) & (
        df['timepoint'] == 'T44')
    df = df[filter]
    fig = px.box(df, x='treatment', y='fixed_total_ratio_number', color='treatment', color_discrete_sequence=colors[abb],
                 points='all', category_orders={'timepoint': ['T11', 'T22', 'T33', 'T44']}, height=h, width=w/3)
    fig.update_traces(boxmean=True, quartilemethod="linear",
                      pointpos=0, jitter=1)
    fig.update_layout(
        showlegend=False, title=strains[s.abbreviations[abb]] + ' transfer 44', title_x=0.5)
    fig.update_xaxes(title='Condition', type='category')
    fig.update_yaxes(title='Proportion of fixed variants', rangemode='tozero')
    fig = font_size(fig, marker_size=5)
    fig.update_yaxes(title_standoff=0)
    fig.update_xaxes(title_standoff=0)
    fig.write_image(join('..', 'plots', 'plots', abb+'_proportion.svg'))


def total_freq_box(y_label, title):
    # Plot for total allele frequency
    hill = get_mutations(add_T0=False)
    # Subsetting for species At and Ct
    a, c = strains[s.abbreviations['at']], strains[s.abbreviations['ct']]
    hill = hill[(hill['strain'] == a) | (hill['strain'] == c)]
    hover_data = ['treatment', 'timepoint', 'strain',  'cosm']

    fig = px.box(hill, x='strain', y='mutations', color='treatment', color_discrete_sequence=colors['all'], facet_col='timepoint',
                 points='all', category_orders={'timepoint': ['T11', 'T22', 'T33', 'T44']},
                 hover_data=hover_data, height=h, width=w)
    fig.update_traces(boxmean=True, quartilemethod="exclusive",
                      pointpos=0, jitter=1)

    # Plot annotations
    titles = ['T11', 'T22', 'T33', 'T44']

    for i, t in enumerate(fig['layout']['annotations']):
        t['text'] = titles[i]

    fig.for_each_yaxis(lambda y: y.update(title=''))
    fig['layout']['yaxis']['title']['text'] = y_label
    fig.for_each_xaxis(lambda y: y.update(title=''))
    fig['layout']['xaxis']['title']['text'] = 'Time-point'
    fig.update_xaxes(title=None)
    fig['layout']['legend']['title']['text'] = 'Condition'
    fig.update_layout(title=title, boxgroupgap=0.2, boxgap=0.3)
    # Setting offsetgroups not ideal
    offsetgroups = ['1', '1', '1', '1',
                    '1', '1', '1', '1',
                    '2', '2', '2', '2',
                    '3', '3', '3', '3']
    for i, d in enumerate(fig['data']):
        d['offsetgroup'] = offsetgroups[i]

    # Setting dticks depending if plotting fixed SNPs or variants

    fig.for_each_yaxis(lambda yaxis: yaxis.update(rangemode="tozero"))

    fig = font_size(fig)
    fig.update_yaxes(title_standoff=0)
    fig.update_xaxes(title_standoff=0)
    fig.write_image(join('..', 'plots', 'plots', 'total_allele_freq_box.svg'))
    return fig


def total_freq_line(abb):
    df = get_mutations(add_T0=False)
    df = df.sort_values(by=['treatment', 'timepoint'], ascending=True)
    df = df[df['strain'] == strains[s.abbreviations[abb]]]
    fig = px.line(df, x='timepoint', y='mutations', width=w/2, height=h,
                  color='treatment', line_group='linegroup', markers=True)
    fig = font_size(fig)
    fig.update_xaxes(title='Transfer')
    fig.update_yaxes(title='Total alllele frequency')
    fig['layout']['legend']['title']['text'] = 'Condition'
    for d in fig['data']:
        d['line']['color'] = colors_t[str(d['name'])]
        d['line']['width'] = 1.5
        d['opacity'] = 0.8
        d['marker']['size'] = 5
    fig.for_each_yaxis(lambda yaxis: yaxis.update(
        tickmode='linear', dtick=2))
    fig.update_xaxes(tickvals=['T11', 'T22', 'T33', 'T44'], ticktext=[
                     '11', '22', '33', '44'])
    fig.update_yaxes(title_standoff=0)
    fig.update_xaxes(title_standoff=0)
    fig.for_each_yaxis(lambda yaxis: yaxis.update(rangemode="tozero"))
    fig.write_image(join('..', 'plots', 'plots',
                    abb+'_total_allele_freq_line.svg'))
    fig = px.line(df, x='timepoint', y='fixed_total_ratio', width=w/2, height=h,
                  color='treatment', line_group='linegroup', markers=True)
    fig = font_size(fig)
    fig.update_xaxes(title='Transfer')
    fig.update_yaxes(title='Ratio fixed variants')
    fig['layout']['legend']['title']['text'] = 'Condition'
    for d in fig['data']:
        d['line']['color'] = colors_t[str(d['name'])]
        d['line']['width'] = 1.5
        d['opacity'] = 0.8
        d['marker']['size'] = 5
    fig.update_xaxes(tickvals=['T11', 'T22', 'T33', 'T44'], ticktext=[
                     '11', '22', '33', '44'])
    fig.update_yaxes(title_standoff=0)
    fig.update_xaxes(title_standoff=0)
    fig.for_each_yaxis(lambda yaxis: yaxis.update(rangemode="tozero"))
    fig.write_image(join('..', 'plots', 'plots',
                    abb+'_fixed_ratio_line.svg'))

    fig = px.line(df, x='timepoint', y='fixed', width=w/2, height=h,
                  color='treatment', line_group='linegroup', markers=True)
    fig = font_size(fig)
    fig.update_xaxes(title='Transfer')
    fig.update_yaxes(title='Number of fixed variants')
    fig['layout']['legend']['title']['text'] = 'Condition'
    for d in fig['data']:
        d['line']['color'] = colors_t[str(d['name'])]
        d['line']['width'] = 1.5
        d['opacity'] = 0.8
        d['marker']['size'] = 5
    fig.update_xaxes(tickvals=['T11', 'T22', 'T33', 'T44'], ticktext=[
                     '11', '22', '33', '44'])
    fig.for_each_yaxis(lambda yaxis: yaxis.update(
        tickmode='linear', dtick=2))
    fig.update_yaxes(title_standoff=0)
    fig.update_xaxes(title_standoff=0)
    fig.for_each_yaxis(lambda yaxis: yaxis.update(rangemode="tozero"))
    ret = fig
    fig.write_image(join('..', 'plots', 'plots',
                    abb+'_fixed_line.svg'))

    f = join('..', 'variants', 'variants_comp_mapping.csv')
    df = get_variants(f, 'illumina')
    df = df[df['strain'] == strains[s.abbreviations[abb]]]
    labels = {'T11': '1', 'T22': '2', 'T33': '3', 'T44': '4'}
    df = df.sort_values(by=['treatment', 'timepoint'], ascending=True)
    fig = px.line(df, x='timepoint', y='hill', width=w/2, height=h, labels=labels,
                  color='treatment', line_group='cosm', markers=True)
    fig = font_size(fig)
    fig.update_xaxes(title='Transfer')
    fig.update_yaxes(title='Number of variants')
    fig['layout']['legend']['title']['text'] = 'Condition'
    for d in fig['data']:
        d['line']['color'] = colors_t[str(d['name'])]
        d['line']['width'] = 1.5
        d['opacity'] = 0.8
        d['marker']['size'] = 5
    fig.update_xaxes(tickvals=['T11', 'T22', 'T33', 'T44'], ticktext=[
                     '11', '22', '33', '44'])
    fig.for_each_yaxis(lambda yaxis: yaxis.update(
        tickmode='linear', dtick=5))
    fig.update_yaxes(title_standoff=0)
    fig.update_xaxes(title_standoff=0)
    fig.for_each_yaxis(lambda yaxis: yaxis.update(rangemode="tozero"))
    fig.write_image(join('..', 'plots', 'plots',
                    abb+'_variants_line.svg'))
    return ret


def fixed_mutations(abb):
    df = get_mutations(add_T0=False)
    df = df.sort_values(by='treatment', ascending=True)
    df = df[df['strain'] == strains[s.abbreviations[abb]]]
    df['treatment'] = df['treatment'].astype(str)
    hover_data = ['cosm', 'treatment', 'timepoint']
    fig = px.scatter(df, x='mutations', y='fixed', width=w, height=h,
                     color='treatment', hover_data=hover_data, color_discrete_sequence=colors[abb], opacity=0.7)
    fig.add_scatter(x=list(range(2, 8)), y=list(range(2, 8)), line={
                    'color': 'gray', 'width': 1}, mode='lines', showlegend=False)
    fig['layout']['legend']['title']['text'] = 'Condition'
    fig.update_xaxes(title='Total alllele frequency')
    fig.update_yaxes(title='Fixed mutations')
    fig.for_each_yaxis(lambda yaxis: yaxis.update(
        tickmode='linear', dtick=2))
    fig.for_each_xaxis(lambda yaxis: yaxis.update(
        tickmode='linear', dtick=2))
    fig = font_size(fig)
    for d in fig['data']:
        d['marker']['size'] = 4
    fig.update_yaxes(title_standoff=0)
    fig.update_yaxes(title_standoff=0)
    fig.update_xaxes(title_standoff=0)
    fig.write_image(join('..', 'plots', 'plots', abb+'_fixed_correlation.svg'))
    return fig


def ratio_fixed_mutations(abb):
    df = get_mutations(add_T0=False)
    df = df.sort_values(by='treatment', ascending=True)
    df = df[(df['strain'] == strains[s.abbreviations[abb]])]
    df['treatment'] = df['treatment'].astype(str)
    hover_data = ['cosm', 'treatment', 'timepoint']

    fig = px.box(df, x='timepoint', y='fixed_total_ratio', height=h, width=w,
                 color='treatment', hover_data=hover_data, color_discrete_sequence=colors[abb],
                 points='all', category_orders={'timepoint': ['T11', 'T22', 'T33', 'T44']})
    fig.update_layout(title='', boxgroupgap=0.2, boxgap=0.3)
    fig.update_traces(boxmean=True, quartilemethod="exclusive",
                      pointpos=0, jitter=1)
    titles = ['T11', 'T22', 'T33', 'T44']

    for i, t in enumerate(fig['layout']['annotations']):
        t['text'] = titles[i]
    offsetgroups = ['1', '1', '1', '1',
                    '1', '1', '1', '1',
                    '2', '2', '2', '2',
                    '3', '3', '3', '3']
    for i, d in enumerate(fig['data']):
        # d['offsetgroup'] = offsetgroups[i]
        pass
    fig['layout']['legend']['title']['text'] = 'Condition'
    fig.update_layout(title='')
    fig.for_each_yaxis(lambda y: y.update(title=''))
    fig['layout']['yaxis']['title']['text'] = 'Ratio fixed mutations'
    fig['layout']['xaxis']['title']['text'] = 'Time-point'
    fig = font_size(fig)
    fig.update_yaxes(title_standoff=0)
    fig.update_xaxes(title_standoff=0)
    fig.write_image(join('..', 'plots', 'plots', 'ratio_fixed_mutations.svg'))
    return fig


def diversity_illumina(f, y_label, title):
    """Plots summed variants or SNPs"""
    hill = get_variants(f, 'illumina')
    # Subsetting for species At and Ct
    a, c = strains[s.abbreviations['at']], strains[s.abbreviations['ct']]
    hill = hill[(hill['strain'] == a) | (hill['strain'] == c)]
    # Color sampling
    hover_data = ['treatment', 'timepoint', 'strain', 'hill', 'cosm']

    fig = px.box(hill, x='strain', y='hill', color='treatment', color_discrete_sequence=colors['all'], facet_col='timepoint',
                 points='all', category_orders={'timepoint': ['T11', 'T22', 'T33', 'T44']},
                 hover_data=hover_data, height=h, width=w)
    fig.update_traces(boxmean=True, quartilemethod="exclusive",
                      pointpos=0, jitter=1)

    # Plot annotations
    titles = ['T11', 'T22', 'T33', 'T44']

    for i, t in enumerate(fig['layout']['annotations']):
        t['text'] = titles[i]

    fig.for_each_yaxis(lambda y: y.update(title=''))
    fig['layout']['yaxis']['title']['text'] = y_label
    fig.update_xaxes(title=None)
    fig['layout']['legend']['title']['text'] = 'Condition'
    fig.update_layout(title=title, boxgroupgap=0.2, boxgap=0.3)
    # Setting offsetgroups not ideal
    offsetgroups = ['1', '1', '1', '1',
                    '1', '1', '1', '1',
                    '2', '2', '2', '2',
                    '3', '3', '3', '3']
    for i, d in enumerate(fig['data']):
        d['offsetgroup'] = offsetgroups[i]

    # Setting dticks depending if plotting fixed SNPs or variants
    if 'snps' in f:
        fig.for_each_yaxis(lambda yaxis: yaxis.update(
            tickmode='linear', dtick=1))
        n = 'number_of_fixed_variants.svg'
    else:
        fig.for_each_yaxis(lambda yaxis: yaxis.update(rangemode="tozero"))
        n = 'number_of_variants.svg'

    fig = font_size(fig)
    fig.update_yaxes(title_standoff=0)
    fig.update_xaxes(title_standoff=0)
    fig.write_image(join('..', 'plots', 'plots',
                    n))
    return fig


def diversity_pacbio():
    """Plots Pacbio SNPs"""
    df = get_variants(
        join('..', 'variants', 'snps_pacbio.csv'), platform='pacbio')

    hover_data = ['hill', 'cosm']
    fig = px.box(df, x='strain', y='hill', color='treatment', color_discrete_sequence=colors['all'], hover_data=hover_data,
                 log_y=False, category_orders={'timepoint': ['T11', 'T22', 'T33', 'T44']}, points='all', width=w, height=h)
    fig['layout']['legend']['title']['text'] = 'Condition'
    fig.update_xaxes(title='Time-point')
    fig.update_yaxes(title='SNPs')
    fig.update_layout(title='')
    fig = font_size(fig)
    fig.update_traces(boxmean=True)
    fig.for_each_yaxis(lambda yaxis: yaxis.update(rangemode="tozero"))
    fig.update_layout(boxgroupgap=0.2, boxgap=0.3)
    fig.update_traces(boxmean=True, quartilemethod="exclusive",
                      pointpos=0, jitter=1)
    offsetgroups = ['1', '1', '2', '3']
    for i, d in enumerate(fig['data']):
        d['offsetgroup'] = offsetgroups[i]
    fig.write_image(join('..', 'plots', 'plots', 'pacbio_variants.svg'))
    return fig


def diversity_ct_box(f, y_label, title):
    """Plots Ct SNPs/variatns with timepoints on x scale"""
    df = get_variants(f, platform='illumina')
    df = df[df['strain'] == strains[s.abbreviations['ct']]]
    fig = px.box(df, x='timepoint', y='hill', color='treatment', color_discrete_sequence=colors['ct'],
                 log_y=False, category_orders={'timepoint': ['T11', 'T22', 'T33', 'T44']}, points='all', width=w, height=h)
    fig['layout']['legend']['title']['text'] = 'Condition'
    fig.update_xaxes(title='Time-point')
    fig.update_yaxes(title=y_label)
    fig.update_layout(title=title)
    fig = font_size(fig)
    fig.update_traces(boxmean=True, quartilemethod="exclusive",
                      pointpos=0, jitter=1)
    if 'snps' in f:
        fig.for_each_yaxis(lambda yaxis: yaxis.update(
            tickmode='linear', dtick=1))
        n = 'number_of_fixed_variants_ct.svg'
    else:
        fig.for_each_yaxis(lambda yaxis: yaxis.update(rangemode="tozero"))
        n = 'number_of_variants_ct.svg'
    fig.write_image(join('..', 'plots', 'plots', n))

    # fig.show()


def mutation_rates(abb, y_label, rates):
    h = 250
    w = 200
    df = get_mutations(add_T0=False)
    df = df[df['strain'] == strains[s.abbreviations[abb]]]
    df = df.sort_values(by='timepoint', ascending=True)
    fig = px.line(df, x='timepoint', y=rates, width=w,
                  height=h, color='treatment', line_group='linegroup', markers=True)
    fig = font_size(fig)
    fig.update_xaxes(title='Time-point')
    fig.update_yaxes(title=y_label)
    fig['layout']['legend']['title']['text'] = 'Condition'
    for d in fig['data']:
        d['line']['color'] = colors_t[str(d['name'])]
        d['line']['width'] = 1.5
        d['opacity'] = 0.7
        d['marker']['size'] = 2
    fig.for_each_yaxis(lambda yaxis: yaxis.update(
        tickmode='linear', dtick=2))
    fig.update_yaxes(title_standoff=0)
    fig.update_xaxes(title_standoff=0)
    f = join('..', 'plots', 'plots', y_label.replace(' ', '_') + '.svg')
    fig.for_each_yaxis(lambda yaxis: yaxis.update(rangemode="tozero"))
    fig.write_image(f)


def trajectories(f, species, title):
    df = pd.read_csv(f, dtype={
                     'cosm': str, 'treatment': str})
    df = df[df['strain'] == s.abbreviations[species]]

    hover_data = ['treatment', 'timepoint', 'cosm', 'depth']
    fig = px.line(df, x='timepoint', y='freq', line_group='linegroup', color_discrete_sequence=colors[species], facet_row='cosm',
                  facet_col='treatment', facet_col_wrap=4, color='treatment', hover_data=hover_data, facet_col_spacing=0.05,
                  category_orders={'timepoint': ['T11', 'T22', 'T33', 'T44']}, markers=True, height=h, width=w)

    conditions = ['Ct condition 2', 'Ct condition 3', 'Ct condition 4']
    cosms = ['M ' + str(i)
             for i in sorted(list(set(df['cosm'])), reverse=True)]
    titles = conditions + cosms
    for i, t in enumerate(fig['layout']['annotations']):
        t['text'] = titles[i]

    fig.update_xaxes(tickvals=['T11', 'T22', 'T33', 'T44'], ticktext=[
                     '11', '22', '33', '44'])
    fig.update_xaxes(title='')
    fig.for_each_yaxis(lambda y: y.update(title=''))
    fig.update_layout(xaxis2=dict(title="Transfer"), yaxis7=dict(
        title='Variant frequency'), showlegend=False)
    fig.update_layout(title=title)
    fig['layout']['legend']['title']['text'] = 'Microcosm'
    fig.for_each_yaxis(lambda yaxis: yaxis.update(
        tickmode='linear', dtick=0.5))
    fig = font_size(fig)
    fig.update_yaxes(title_standoff=0)
    fig.update_xaxes(title_standoff=0)
    fig.write_image(join('..', 'plots', 'plots',
                    species+'_trajectories.svg'))
    return fig


def t_test(df, column):
    """T test for testing variant richness significance."""
    at = df[df['strain'] == strains[s.abbreviations['at']]]
    for j in ['T11', 'T22', 'T33', 'T44']:
        g1 = at[(at['timepoint'] == j) & (
            at['treatment'] == 3)][column].to_numpy()
        g2 = at[(at['timepoint'] == j) & (
            at['treatment'] == 4)][column].to_numpy()
        t, p = ttest_ind(g1, g2)
        if p < 0.05:
            print('At', j)
            print('T-Statistic', t, 'P-Value', p)

    ct = df[df['strain'] == strains[s.abbreviations['ct']]]
    for j in ['T11', 'T22', 'T33', 'T44']:
        g1 = ct[(ct['timepoint'] == j) & (
            ct['treatment'] == 2)][column].to_numpy()
        g2 = ct[(ct['timepoint'] == j) & (
            ct['treatment'] == 3)][column].to_numpy()
        g3 = ct[(ct['timepoint'] == j) & (
            ct['treatment'] == 4)][column].to_numpy()
        t, p = kruskal(g1, g2)
        if p < 0.05:
            print('Ct treatment 2 and 3', j)
            print('T-Statistic', t, 'P-Value', p)
        t, p = kruskal(g1, g3)
        if p < 0.05:
            print('Ct treatment 2 and 4', j)
            print('T-Statistic', t, 'P-Value', p)
        t, p = kruskal(g2, g3)
        if p < 0.05:
            print('Ct treatment 3 and 4', j)
            print('T-Statistic', t, 'P-Value', p)


def plotter():
    variants = join('..', 'variants', 'variants_comp_mapping.csv')
    # total_freq_line('at')
    fig = trajectories(variants, 'ct', '')
    variants = join('..', 'variants', 'variants_comp_mapping.csv')
    snps = join('..', 'variants', 'snps_freebayes_comp_mapping.csv')
    total_freq_box('Total allele frequency', '')
    total_freq_line('ct')
    fixed_mutations('ct')
    ratio_fixed_mutations('ct')
    variants = join('..', 'variants', 'variants_comp_mapping.csv')
    fig = diversity_illumina(variants, 'Number of variants', '')
    fig = diversity_illumina(snps, 'Number of fixed mutations', '')
    diversity_ct_box(variants, 'Number of variants', '')
    trajectories(variants, 'ct', '')
    mutation_rates('ct', 'Fixation rate', 'fixation_rate')
    mutation_rates('ct', 'Accumulation rate', 'acc_rate')
    coverage()


def generation_time():
    out = pd.DataFrame(columns=['treatment', 'microcosm',
                                'generation_time', 'time', 'lg'])
    df = pd.read_csv(join('..', 'variants', 'cfus_ct.csv'))
    for c in df.columns:
        for i in df.index[:-1]:
            try:
                Nt = df[c][i + 1] * 100
                N0 = df[c][i]
                gt = 7*24 / (math.log2(Nt / N0))
            except ValueError:
                gt = None
            except ZeroDivisionError:
                gt = None
            t, m = c.split('.')[1], c.split('.')[2]
            out.loc[len(out)] = [t, m, gt, i, c]
    out = out[(out['treatment'] == '2') | (out['treatment'] == '4')]
    colors = px.colors.sample_colorscale(
        "Agsunset", [n/(2 - 1) for n in range(2)])
    fig = px.box(out, x='treatment', y='generation_time', points='all', color_discrete_sequence=colors,
                 color='treatment')
    fig.update_traces(boxmean=True, quartilemethod="exclusive",
                      pointpos=0, jitter=1)
    fig.show()


def ct_growth_curves():
    df = pd.read_csv('../variants/cfus_ct.csv')
    tmp = pd.DataFrame(columns=['Condition', 'lg', 'time', 'count'])

    for c in df.columns:
        for i, j in enumerate(df[c]):
            tmp.loc[len(tmp)] = [c[3], c, i, j]
    tmp = tmp.astype({'count': int})
    fig = px.line(tmp, x='time', y='count', line_group='lg',
                  color='Condition', log_y=True, width=w, height=h)
    for d in fig['data']:
        d['line']['color'] = colors_t[d['name']]
    fig.update_xaxes(title='Transfer')
    fig.update_yaxes(title='CFUs/mL')
    fig.update_layout(yaxis=dict(exponentformat="E"))
    fig = font_size(fig, line_size=1)
    fig.write_image(join('..', 'plots', 'plots', 'ct_growth_curves.svg'))


def annot_chrom(abb, chromosome, length, width, fname, cosm_names, fsize=10, height=h):
    df = pd.read_csv(join('..', 'annotations', abb +
                     '_variants_annotations.csv'))
    mask = (df['strain'] == s.abbreviations[abb]) & (
        df['chrom'] == chromosome) & (df['timepoint'] == 'T44')
    df = df[mask]
    df = df.astype({'treatment': str, 'cosm': str})
    df['freq_color'] = df['freq']
    df['freq'] = 0
    df['facet_row'] = df['treatment'] + '_' + df['cosm']
    df.insert(len(df.columns), 'coding', None)
    for i, row in df.iterrows():
        if row['gene'] == 'Not annotated':
            df.at[i, 'coding'] = False
        else:
            df.at[i, 'coding'] = True
    df = df.sort_values(by='facet_row')
    lines = {'at': [(1, 1), (1, 2), (3, 1), (3, 2), (3, 3), (3, 4), (3, 5), (4, 1), (4, 2), (4, 3), (4, 4), (4, 5)],
             'ct': [(2, 1), (2, 2), (2, 3), (2, 4), (2, 5), (3, 1), (3, 2), (3, 3), (3, 4), (3, 5), (4, 1), (4, 2), (4, 3), (4, 4), (4, 5)]}
    row_titles = []
    if cosm_names:
        for t, m in lines[abb]:
            row_titles.append('M'+str(m))
    fig = make_subplots(
        rows=len(lines[abb]), cols=1, shared_xaxes=True, row_titles=row_titles)
    for i, j in enumerate(lines[abb]):
        fig.add_trace(go.Scatter(x=[0, length], y=[
                      0, 0], mode='lines', line=dict(color='gray')), col=1, row=i+1)

    for i, (t, m) in enumerate(lines[abb]):
        mask = (df['treatment'] == str(t)) & (df['cosm'] == str(m))
        tmp = df[mask]
        for j, row in tmp.iterrows():
            if row['coding']:
                symbol = 'x-thin'
            else:
                symbol = 'circle-open'

            x, y = [row['pos']], [row['freq']]
            fig.add_trace(go.Scatter(x=x, y=y, name=str(t)+'_'+str(m), mode='markers', marker=dict(color='black',
                                                                                                   symbol=symbol, line=dict(width=0.8, color='black'))
                                     ), row=i+1, col=1)

    fig.for_each_yaxis(lambda axis: axis.update(visible=False))
    fig.for_each_xaxis(lambda axis: axis.update(showgrid=False))
    fig.for_each_xaxis(lambda axis: axis.update(zeroline=False))
    fig.for_each_annotation(lambda axis: axis.update(textangle=0))
    fig.update_xaxes(rangemode='tozero')
    symbols = {'True': 'x-thin',
               'False': 'circle-open'}

    fig = font_size(fig, marker_size=3, line_size=4, fsize=fsize)

    fig.update_layout(
        plot_bgcolor='rgba(0,0,0,0)',  # Fully transparent background
        paper_bgcolor='rgba(0,0,0,0)',
        showlegend=False,
        autosize=False,
        minreducedwidth=width,
        minreducedheight=height,
        width=width,
        height=height,
    )
    # fig.data = fig.data[::-1]

    if abb == 'ct':
        for d in fig['data'][:5]:
            d['line']['color'] = colors_t['2']
        for d in fig['data'][5:10]:
            d['line']['color'] = colors_t['3']
        for d in fig['data'][10:15]:
            d['line']['color'] = colors_t['4']
    if abb == 'at':
        for d in fig['data'][:2]:
            d['line']['color'] = colors_t['1']
        for d in fig['data'][2:7]:
            d['line']['color'] = colors_t['3']
        for d in fig['data'][7:12]:
            d['line']['color'] = colors_t['4']
    fig.write_image(join('..', 'plots', 'plots', fname))
    return fig


def annotater():
    fig = annot_chrom('ct', 'ct_0', 5931804, w/5*4,
                      'ct_chrom_annot.svg', False)
    fig = annot_chrom('ct', 'ct_1', 198853, w/5, 'ct_plas_annot.svg', True)
    fig = annot_chrom('at', 'at_0', 3000155, w/7*2,
                      'at_chrom_annot.svg', False, fsize=5, height=h+3.5)
    fig = annot_chrom('at', 'at_1', 1955452, w/7*2,
                      'at_plas1_annot.svg', False, fsize=5, height=h+3.5)
    fig = annot_chrom('at', 'at_2', 214247, w/7,
                      'at_plas2_annot.svg', False, fsize=5, height=h+3.5)
    fig = annot_chrom('at', 'at_3', 229513, w/7,
                      'at_plas3_annot.svg', False, fsize=6, height=h+3.5)
    fig = annot_chrom('at', 'at_4', 31100, w/7,
                      'at_plas4_annot.svg', True, fsize=6, height=h+3.5)


def zero_coverage():
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
    fig.write_image('tmp.svg')
