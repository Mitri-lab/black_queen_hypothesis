from scipy.stats import ttest_ind, kruskal
from samples import Samples
from os.path import join, exists
import vcfpy
import pandas as pd
import plotly.express as px
from plotly.subplots import make_subplots
import math
from scipy.spatial import distance
# Sample class for parsing contact me for infos
# eric.ulrich@unil.ch
s = Samples()

# Global names for species used for plotting
strains = {s.abbreviations['at']: 'At',
           s.abbreviations['ct']: 'Ct',
           s.abbreviations['oa']: 'Oa',
           s.abbreviations['ms']: 'Ms'}


colors = {'ct': ['#a431a0', '#ea4f88', '#f89178'],
          'at': ['#4b2991', '#ea4f88', '#f89178'],
          'all': ['#4b2991', '#a431a0', '#ea4f88', '#f89178'],
          'cosm': ['#4b2991', '#a431a0', '#ea4f88', '#f89178', '#edd9a3']
          }

colors_t = {'1': '#4b2991',
            '2': '#a431a0',
            '3': '#ea4f88',
            '4': '#f89178'}


# This section is executed on the cluster for parsing data on cluster
w = 400
h = 250

def parse_variants():
    """Parses vcf files and dumps variants with metdadata as csv."""
    dfs = []
    for strain, samples in s.strains.items():
        for sample in samples:
            if sample['platform'] == 'illumina':
                f = join(sample['dir_name'], 'var.annot.vcf')
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
                            snp['total_freq'] = sum(
                                record.INFO['AO']) / record.INFO['DP']
                            snp['alt'] = r
                            snp['alt_count'] = record.INFO['AO'][i]
                            snp['ref'] = record.REF
                            snp['type'] = record.INFO['TYPE'][i]
                            snp['len'] = record.INFO['LEN'][i]
                            snp['eff'] = record.INFO['ANN'][i].split('|')[1]
                            # Creates unique key per variant used for plotting trajectories
                            key = '.'.join([snp['chrom'], str(snp['pos']),
                                            str(sample['treatment']), str(sample['cosm']), str(snp['alt'])])
                            snp['linegroup'] = key
                            # Quality cutoff
                            if snp['qual'] >= 20:
                                if (sample['strain'] == s.abbreviations['oa']) & (snp['chrom'] == 'tig00000002_polypolish') & (snp['pos'] in range(119960, 120700)):
                                    pass
                                else:
                                    df = (pd.concat(
                                        [pd.DataFrame(snp, index=[0]), pd.DataFrame(sample, index=[0])], axis=1))
                                    dfs.append(df)
    out = pd.concat(dfs)
    out['treatment'] = out['treatment'].astype(str)
    out['cosm'] = out['cosm'].astype(str)
    out.to_csv(join('..', 'variants', 'variants_comp_mapping.csv'), index=False)
    filter = (out['alt_count'] >= 3) & (out['freq'] >= 0.1) & (
        out['qual'] >= 20) & (out['eff'] != 'synonymous_variant')
    out[filter].to_csv(
        join('..', 'variants', 'ns_variants_comp_mapping.csv'), index=False)
    filter = (out['alt_count'] >= 3) & (out['freq'] >= 0.1) & (
        out['qual'] >= 20)
    out[filter].to_csv(
        join('..', 'variants', 'variants_comp_mapping.csv'), index=False)
    filter = (out['alt_count'] >= 3) & (out['total_freq'] >= 0.95) & (
        out['qual'] >= 20) & (out['eff'] != 'synonymous_variant')
    out[filter].to_csv(
        join('..', 'variants', 'ns_snps_freebayes_comp_mapping.csv'), index=False)
    filter = (out['alt_count'] >= 3) & (out['total_freq'] >= 0.95) & (
        out['qual'] >= 20)
    out[filter].to_csv(
        join('..', 'variants', 'snps_freebayes_comp_mapping.csv'), index=False)


def parse_snps():
    """Parses fixed SNPs from snippy"""
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
    out.to_csv(join('..', 'variants', 'snps_comp_mapping.csv'), index=False)
    out = out[out['eff'] != 'synonymous_variant']
    out.to_csv(join('..', 'variants', 'ns_snps_comp_mapping.csv'), index=False)


def parse_snps_pacbio():
    """Parses SNPs from pacbio data.
    Clonal data, SNPs identified with assemblies.
    Data equal to strain data."""
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
    """Parses coverage textfiles created with snakemake workflow."""
    df = pd.DataFrame(columns=['sample', 'timepoint',
                               'species', 'treatment', 'depth', 'cosm'])
    i = 0
    for strain, samples in s.strains.items():
        for sample in samples:
            if sample['platform'] == 'illumina':
                f = join(sample['dir_name'], 'coverage.txt')
                with open(f, 'r') as handle:
                    d = float(handle.read())
                    if d < 1:
                        d = 0
                r = [sample['name'], sample['timepoint'], strains[strain],
                     sample['treatment'], d, sample['cosm']]
                df.loc[i] = r
                i += 1
    df.to_csv(join('..', 'variants', 'mean_coverage.csv'), index=False)

# Local execution starts here for plotting


def font_size(fig):
    """Style function for figures setting fot size and true black color."""
    for d in fig['data']:
        d['marker']['size'] = 3
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


def get_variants(f, platform):
    """Sums up SNPs or variants from variant csv files in ../variants"""
    snps = pd.read_csv(f)
    # Hill name from hill numbers, q = 0. Artifiact from old analysis still valid.
    # snps = snps[snps['eff'] != 'synonymous_variant']
    hill = pd.DataFrame(
        columns=['strain', 'treatment', 'hill', 'timepoint', 'cosm'])
    i = 0
    for strain, samples in s.strains.items():
        for sample in samples:
            if sample['platform'] == platform:
                # Subsetting dataframe per species per sample
                sub = snps[(snps['strain'] == strain) & (
                    snps['name'] == sample['name'])]
                hill.loc[i] = [strains[strain], sample['treatment'],
                               len(sub), sample['timepoint'], sample['cosm']]

                i += 1
    hill = hill.sort_values(by='treatment', ascending=True)
    return hill


def get_mutations(add_T0=True):
    variants = pd.read_csv(
        join('..', 'variants', 'variants_comp_mapping.csv'))
    snps = pd.read_csv(
        join('..', 'variants', 'snps_freebayes_comp_mapping.csv'))
    out = pd.DataFrame(columns=['strain', 'name', 'cosm',
                                'treatment', 'timepoint', 'mutations', 'fixed', 'linegroup'])
    for strain, samples in s.strains.items():
        for sample in samples:
            if sample['platform'] == 'illumina':
                lg = '_'.join([s.abbreviations[strain], str(
                    sample['cosm']), str(sample['treatment'])])
                tmp_var = variants[(variants['strain'] == strain) & (
                    variants['name'] == sample['name'])]
                fixed = snps[(snps['strain'] == strain) & (
                    snps['name'] == sample['name'])]
                out.loc[len(out)] = [strains[strain], sample['name'], sample['cosm'],
                                     sample['treatment'], sample['timepoint'], sum(tmp_var['freq']), len(fixed),  lg]
                if add_T0:
                    out.loc[len(out)] = [strains[strain], sample['name'],
                                         sample['cosm'], sample['treatment'], 'T0', 0, 0, lg]
    for i, (f, m) in enumerate(zip(out['fixed'], out['mutations'])):
        try:
            out.at[i, 'fixed_total_ratio'] = f/m
        except ZeroDivisionError:
            out.at[i, 'fixed_total_ratio'] = None
    out.insert(len(out.columns), 'fixation_rate', None)
    out.insert(len(out.columns), 'acc_rate', None)

    for i, row in out.iterrows():
        t, fixed, mutations = row['timepoint'], row['fixed'], row['mutations']
        t = int(t[1])
        out.at[i, 'fixation_rate'] = fixed / t
        out.at[i, 'acc_rate'] = mutations / t
    out = out.sort_values(by='treatment', ascending=True)
    out['hill'] = out['mutations']
    out.to_csv(join('..', 'variants', 'total_allele_frequncies.csv'), index=False)
    return out


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
    df = df.sort_values(by=['treatment','timepoint'], ascending=True)
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
        d['x'] = [1,2,3,4]
    fig.for_each_yaxis(lambda yaxis: yaxis.update(
        tickmode='linear', dtick=2))
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
        d['x'] = [1,2,3,4]
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
        d['x'] = [1,2,3,4]
    fig.for_each_yaxis(lambda yaxis: yaxis.update(
        tickmode='linear', dtick=2))
    fig.update_yaxes(title_standoff=0)
    fig.update_xaxes(title_standoff=0)
    fig.for_each_yaxis(lambda yaxis: yaxis.update(rangemode="tozero"))
    fig.write_image(join('..', 'plots', 'plots',
                    abb+'_fixed_line.svg'))

    f = join('..', 'variants', 'variants_comp_mapping.csv')
    df = get_variants(f, 'illumina')
    df = df[df['strain'] == strains[s.abbreviations[abb]]]
    labels = {'T11':'1','T22':'2','T33':'3','T44':'4'}
    df = df.sort_values(by=['treatment','timepoint'], ascending=True)
    fig = px.line(df, x='timepoint', y='hill', width=w/2, height=h,labels=labels,
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
        d['x'] = [1,2,3,4]
    fig.for_each_yaxis(lambda yaxis: yaxis.update(
        tickmode='linear', dtick=5))
    fig.update_yaxes(title_standoff=0)
    fig.update_xaxes(title_standoff=0)
    fig.for_each_yaxis(lambda yaxis: yaxis.update(rangemode="tozero"))
    fig.write_image(join('..', 'plots', 'plots',
                    abb+'_variants_line.svg'))
    return fig


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


def coverage():
    df = pd.read_csv(join('..', 'variants', 'mean_coverage.csv'))
    # Possibility to subset for all species
    # a, c = strains[s.abbreviations['at']], strains[s.abbreviations['ct']]
    # df = df[(df['species'] == a) | (df['species'] == c)]
    df = df.sort_values(by='treatment', ascending=True)
    depth = [None if i == 0 else i for i in df['depth']]
    df['depth'] = depth
    hover_data = ['treatment', 'timepoint', 'species', 'depth', 'cosm']
    fig = px.box(df, x='species', y='depth', color='treatment', facet_col='timepoint', points='all',
                 category_orders={'timepoint': [
                       'T11', 'T22', 'T33', 'T44']}, color_discrete_sequence=colors['all'],
                 log_y=True, hover_data=hover_data, height=h, width=w)
    titles = ['T11', 'T22', 'T33', 'T44']
    for i, t in enumerate(fig['layout']['annotations']):
        t['text'] = titles[i]

    fig.for_each_yaxis(lambda y: y.update(title=''))
    fig['layout']['yaxis']['title']['text'] = 'Mean coverage'
    fig.update_xaxes(title=None)
    fig['layout']['legend']['title']['text'] = 'Condition'
    fig.update_layout(title='', boxgroupgap=0.2, boxgap=0.3)
    fig.update_traces(boxmean=True, quartilemethod="exclusive",
                      pointpos=0, jitter=1)
    #fig.for_each_yaxis(lambda yaxis: yaxis.update(rangemode="tozero"))

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
                    'coverage_2.svg'))
    return fig
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
    fig = px.line(df, x='timepoint', y='freq', line_group='linegroup', color_discrete_sequence=colors['cosm'],
                  facet_col='treatment', facet_col_wrap=4, color='cosm', hover_data=hover_data, facet_col_spacing=0.05,
                  category_orders={'timepoint': ['T11', 'T22', 'T33', 'T44']}, markers=True, height=h, width=w)

    fig.for_each_yaxis(lambda y: y.update(title=''))
    fig['layout']['yaxis']['title']['text'] = 'Variant frequency'
    titles = ['Condition ' + str(i)
              for i in sorted(list(set(df['treatment'])))]
    for i, t in enumerate(fig['layout']['annotations']):
        t['text'] = titles[i]
    fig.update_xaxes(title='')
    fig.update_layout(xaxis2=dict(title="Transfer"))
    fig.update_layout(title=title)
    fig['layout']['legend']['title']['text'] = 'Microcosm'
    fig.for_each_yaxis(lambda yaxis: yaxis.update(
        tickmode='linear', dtick=0.2))
    fig = font_size(fig)
    fig.update_yaxes(title_standoff=0)
    fig.update_xaxes(title_standoff=0)
    for d in fig['data']:
        d['x'] = [1,2,3,4]
    fig.write_image(join('..', 'plots', 'plots',
                    species+'_trajectories.svg'))


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


variants = join('..', 'variants', 'variants_comp_mapping.csv')
#total_freq_line('ct')

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

f = join('..', 'variants', 'variants_comp_mapping.csv')
df = pd.read_csv(f)
filter = (df['strain'] == s.abbreviations['ct']) & (df['timepoint'] == 'T44')
df = df[filter]
mask = []
for alt in df['alt']:
    if 'SNV' in alt:
        mask.append(True)
    else:
        mask.append(False)
df = df[mask]

vars = []
lgs = list(set(df['linegroup']))
for v in lgs:
    if max(df[df['linegroup'] == v]['freq']) >= 0.2:
        vars.append(v)
d = [[0] * len(vars)]
sample_name = ['Ancestor']
for sample in s.strains[s.abbreviations['ct']]:
    if (sample['platform'] == 'illumina') & (sample['timepoint'] == 'T44'):    
        tmp = df[df['name'] == sample['name']]
        row = []
        for v in vars:
            i = tmp[tmp['linegroup'] == v]
            if len(i) == 0:
                row.append(0)
            elif len(i) == 1:
                row.append(i['freq'])
            else:
                print(i)
        d.append(row)
        sample_name.append(sample['name'])
import numpy as np
from scipy.spatial import distance
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

d = np.array(d,dtype='float64')
distance_matrix = distance.cdist(d, d, 'hamming')
df = pd.DataFrame(distance_matrix)
#df = df.replace(0,3)
df.index,df.columns = sample_name,sample_name
fig = px.imshow(df)
fig.show()
out = pd.DataFrame(columns=['treatment','distance'])
row = df.loc['Ancestor']
row.index = sample_name
for i,j in row.items():
    t = i[4]
    if t == 's':
        pass
    else:
        out.loc[len(out)] = [t,j]

fig = px.box(out,x='treatment',y='distance',points='all')
