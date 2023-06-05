from scipy.stats import ttest_ind
from samples import Samples
from os.path import join, exists
import vcfpy
import pandas as pd
import plotly.express as px
from plotly.subplots import make_subplots

# Sample class for parsing contact me for infos
# eric.ulrich@unil.ch
s = Samples()

# Global names for species used for plotting
strains = {s.abbreviations['at']: 'At',
           s.abbreviations['ct']: 'Ct',
           s.abbreviations['oa']: 'Oa',
           s.abbreviations['ms']: 'Ms'}

# This section is executed on the cluster for parsing data on cluster


def parse_variants():
    """Parses vcf files and dumps variants with metdadata as csv."""
    dfs = []
    for strain, samples in s.strains.items():
        for sample in samples:
            if sample['platform'] == 'illumina':
                f = join(sample['dir_name'], 'var.vcf')
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
    out = out[out['freq'] != 0]
    out['treatment'] = out['treatment'].astype(str)
    out['cosm'] = out['cosm'].astype(str)
    out.to_csv(join('..', 'variants_comp_mapping.csv'), index=False)


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
    out.to_csv(join('..', 'variants', 'snps_comp_mapping.csv'), index=False)


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
    out.to_csv(join('..', 'snps_pacbio.csv'), index=False)


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
    return fig


def get_variants(f, platform):
    """Sums up SNPs or variants from variant csv files in ../variants"""
    snps = pd.read_csv(f)
    # Hill name from hill numbers, q = 0. Artifiact from old analysis still valid.
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


def diversity_illumina(f, y_label, title):
    """Plots summed variants or SNPs"""
    hill = get_variants(f, 'illumina')
    # Subsetting for species At and Ct
    a, c = strains[s.abbreviations['at']], strains[s.abbreviations['ct']]
    hill = hill[(hill['strain'] == a) | (hill['strain'] == c)]
    # Color sampling
    n_colors = len(set(hill['treatment']))
    colors = px.colors.sample_colorscale(
        "Agsunset", [n/(n_colors - 1) for n in range(n_colors)])
    hover_data = ['treatment', 'timepoint', 'strain', 'hill', 'cosm']

    fig = px.box(hill, x='strain', y='hill', color='treatment', color_discrete_sequence=colors, facet_col='timepoint',
                 points='all', category_orders={'timepoint': ['T11', 'T22', 'T33', 'T44']},
                 hover_data=hover_data, height=250, width=400)
    fig.update_traces(boxmean=True, quartilemethod="exclusive",
                      pointpos=0, jitter=1)

    # Plot annotations
    titles = ['T11', 'T22', 'T33', 'T44']

    for i, t in enumerate(fig['layout']['annotations']):
        t['text'] = titles[i]

    fig.for_each_yaxis(lambda y: y.update(title=''))
    fig['layout']['yaxis']['title']['text'] = y_label
    fig.update_xaxes(title=None)
    fig['layout']['legend']['title']['text'] = 'Treatment'
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
        n = 'fixed_variants.svg'
    else:
        fig.for_each_yaxis(lambda yaxis: yaxis.update(rangemode="tozero"))
        n = 'variants.svg'

    fig = font_size(fig)
    fig.write_image(join('..', 'plots', 'snps_figures',
                    n))
    return fig


def ct_box(f, y_label, title):
    """Plots Ct SNPs/variatns with timepoints on x scale"""
    df = get_variants(
        join('..', 'variants', 'variants_comp_mapping.csv'), platform='illumina')
    df = df[df['strain'] == strains[s.abbreviations['ct']]]
    n_colors = len(set(df['treatment']))
    colors = px.colors.sample_colorscale(
        "Agsunset", [n/(n_colors - 1) for n in range(n_colors)])
    fig = px.box(df, x='timepoint', y='hill', color='treatment', color_discrete_sequence=colors,
                 log_y=False, category_orders={'timepoint': ['T11', 'T22', 'T33', 'T44']}, points='all', width=350, height=300)
    fig['layout']['legend']['title']['text'] = 'Treatment'
    fig.update_xaxes(title='Timepoint')
    fig.update_yaxes(title=y_label)
    fig.update_layout(title=title)
    fig = font_size(fig)
    fig.update_traces(boxmean=True)
    if 'snps' in f:
        fig.for_each_yaxis(lambda yaxis: yaxis.update(
            tickmode='linear', dtick=1))
        n = 'fixed_variants_ct.svg'
    else:
        fig.for_each_yaxis(lambda yaxis: yaxis.update(rangemode="tozero"))
        n = 'variants_ct.svg'
    fig.write_image(join('..', 'plots', 'snps_figures',n))


def snp_distribution(f, abb, timepoint, subset=False, add_clusters=False):
    """Plots variant distributions over genome.
    Allows to annotate "clusters". Subsetting possible with treatment
    as string."""
    hill = pd.read_csv(f, dtype={'cosm': str, 'treatment': str})
    hill = hill[hill['strain'] == s.abbreviations[abb]]
    hill = hill[hill['timepoint'] == timepoint]
    hill.index = range(len(hill))
    hill.insert(0, 'ID', hill.index)
    if subset:
        hill = hill[hill['treatment'] == subset]

    hover_data = ['pos', 'depth', 'qual', 'cosm', 'ID']

    if subset:
        n_colors = len(set(hill['cosm']))
        colors = px.colors.sample_colorscale(
            "Agsunset", [n/(n_colors - 1) for n in range(n_colors)])
        fig = px.scatter(hill, x='pos', y='freq', color='cosm', color_discrete_sequence=colors,
                         facet_col='chrom', facet_col_wrap=1,
                         facet_row_spacing=0.12, hover_data=hover_data, width=400, height=300)
    else:
        n_colors = len(set(hill['cosm']))
        colors = px.colors.sample_colorscale(
            "Agsunset", [n/(n_colors - 1) for n in range(n_colors)])
        fig = px.scatter(hill, x='pos', y='freq', color='treatment', color_discrete_sequence=colors,
                         facet_col='chrom', facet_col_wrap=1,
                         facet_row_spacing=0.12, hover_data=hover_data, width=400, height=300)
    fig.update_xaxes(matches=None, showticklabels=True)
    fig.for_each_yaxis(lambda y: y.update(title=''))
    fig['layout']['yaxis']['title']['text'] = 'Variant frequency'
    fig['layout']['xaxis']['title']['text'] = 'Position'
    for i, t in enumerate(fig['layout']['annotations']):
        t['text'] = t['text'].replace('chrom=', '')

    if add_clusters:
        # Annotations for plot
        annot = pd.read_csv('clusters_'+abb+'.csv')
        for i, chrom in enumerate(sorted(set(annot['chrom']))):
            annot_sub = annot[annot['chrom'] == chrom]
            for id in set(annot_sub['id']):
                asid = annot_sub[annot_sub['id'] == id]
                asid.index = range(len(asid))
                x0 = asid.loc[0].pos
                x1 = asid.loc[len(asid)-1].pos
                fig.add_vline(x=(x0+x1)/2, row=i, annotation_text=id,
                              opacity=0.25)
    if subset:
        n = abb+'_positions_cosm.svg'
        fig.update_layout(
            title='Variants colored by microcosm for treatment 4 in '+strains[s.abbreviations[abb]])
        fig['layout']['legend']['title']['text'] = 'Microcosm'

    else:
        fig.update_layout(title='Variant distribution along the genome')
        n = abb+'_positions.svg'
        fig['layout']['legend']['title']['text'] = 'Treatment'
    fig = font_size(fig)
    fig.for_each_yaxis(lambda yaxis: yaxis.update(
        tickmode='linear', dtick=0.2))
    fig.write_image(join('..', 'plots', 'snps_figures',
                    n))
    return fig


def coverage():
    df = pd.read_csv(join('..', 'variants', 'mean_coverage.csv'))
    # Possibility to subset for all species
    #a, c = strains[s.abbreviations['at']], strains[s.abbreviations['ct']]
    #df = df[(df['species'] == a) | (df['species'] == c)]
    df = df.sort_values(by='treatment', ascending=True)
    hover_data = ['treatment', 'timepoint', 'species', 'depth', 'cosm']
    n_colors = len(set(df['treatment']))
    colors = px.colors.sample_colorscale(
        "Agsunset", [n/(n_colors - 1) for n in range(n_colors)])
    fig = px.box(df, x='species', y='depth', color='treatment', facet_col='timepoint', points='all',
                 category_orders={'timepoint': [
                       'T11', 'T22', 'T33', 'T44']}, color_discrete_sequence=colors,
                 log_y=True, hover_data=hover_data, height=250, width=450)
    titles = ['T11', 'T22', 'T33', 'T44']
    for i, t in enumerate(fig['layout']['annotations']):
        t['text'] = titles[i]

    fig.for_each_yaxis(lambda y: y.update(title=''))
    fig['layout']['yaxis']['title']['text'] = 'Mean coverage'
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
    fig.write_image(join('..', 'plots', 'snps_figures',
                    'coverage.svg'))
    return fig
    # fig.show()


def trajectories(species, title):
    df = pd.read_csv(join('..', 'variants', 'variants_comp_mapping.csv'), dtype={
                     'cosm': str, 'treatment': str})
    df = df[df['strain'] == s.abbreviations[species]]
    n_colors = len(set(df['cosm']))
    colors = px.colors.sample_colorscale(
        "Agsunset", [n/(n_colors - 1) for n in range(n_colors)])
    hover_data = ['treatment', 'timepoint', 'cosm', 'depth']
    fig = px.line(df, x='timepoint', y='freq', line_group='linegroup', color_discrete_sequence=colors,
                  facet_col='treatment', facet_col_wrap=4, color='cosm', hover_data=hover_data, facet_col_spacing=0.05,
                  category_orders={'timepoint': ['T11', 'T22', 'T33', 'T44']}, markers=True, height=250, width=400)

    fig.for_each_yaxis(lambda y: y.update(title=''))
    fig['layout']['yaxis']['title']['text'] = 'Variant frequency'
    titles = ['Treatment ' + str(i)
              for i in sorted(list(set(df['treatment'])))]
    for i, t in enumerate(fig['layout']['annotations']):
        t['text'] = titles[i]
    fig.update_xaxes(title='')
    fig.update_layout(xaxis2=dict(title="Timepoint"))
    fig.update_layout(title=title)
    fig['layout']['legend']['title']['text'] = 'Microcosm'
    fig.for_each_yaxis(lambda yaxis: yaxis.update(
        tickmode='linear', dtick=0.2))
    fig = font_size(fig)
    fig.write_image(join('..', 'plots', 'snps_figures',
                    species+'_trajectories.svg'))
    return fig


def cluster_trajectories(abb):
    """This allows to follow the trajectories of SNP clusters of interest."""
    df = pd.read_csv(join('..', 'variants', 'variants_comp_mapping.csv'), dtype={
                     'cosm': str, 'treatment': str})
    df = df[df['strain'] == s.abbreviations[abb]]
    cluster = pd.read_csv('clusters_'+abb+'.csv')
    out = pd.DataFrame(columns=['treatment', 'chrom', 'pos', 'cosm',
                                'timepoint', 'freq', 'qual', 'linegroup', 'cluster', 'strain'])
    k = 0
    for i, c, p in zip(cluster['id'], cluster['chrom'], cluster['pos']):
        df_c = df[df['chrom'] == c]
        df_p = df_c[df_c['pos'] == p]
        for j, row in df_p.iterrows():
            out.loc[k] = [row['treatment'], row['chrom'], row['pos'],
                          row['cosm'], row['timepoint'], row['freq'], row['qual'], row['linegroup'], i, 'ct']
            k += 1
    n_colors = len(set(cluster['id']))
    if n_colors == 1:
        colors = px.colors.sample_colorscale(
            "Agsunset", [1])
    else:
        colors = px.colors.sample_colorscale(
            "Agsunset", [n/(n_colors - 1) for n in range(n_colors)])
    fig = px.line(out, x='timepoint', y='freq', line_group='linegroup',
                  facet_col='treatment', facet_col_wrap=4, color='cluster', color_discrete_sequence=colors, width=400, height=350,
                  category_orders={'timepoint': ['T11', 'T22', 'T33', 'T44'], 'treatment': sorted(set(out['treatment']))}, markers=True)
    titles = ['Treatment ' + str(i)
              for i in sorted(list(set(df['treatment'])))]
    for i, t in enumerate(fig['layout']['annotations']):
        t['text'] = titles[i]
    fig.for_each_yaxis(lambda y: y.update(title=''))
    fig['layout']['yaxis']['title']['text'] = 'Variant frequency'
    fig.update_xaxes(title=None)
    fig.for_each_xaxis(lambda x: x.update(title=''))
    fig['layout']['xaxis']['title']['text'] = 'Variant frequency'
    fig = font_size(fig)
    fig.show()


def annotate_clusters(abb):
    """Prints annotations for identified SNPs."""
    def annotate_pos(gbk, c, p):
        for (start, end), (gene, product) in gbk[c].items():
            if p in range(start, end):
                return [gene, product]
        return False

    f = join('..', 'annotations', abb + '.tsv')
    df = pd.read_csv(f, sep='\t')
    # For plotting we renamed contigs to at_0 etc.
    # Rename of contigs in annotations for hashing.
    contigs = {c: abb+'_'+str(i)
               for i, c in enumerate(sorted(set(df['Sequence Id'])))}
    for i, chrom in enumerate(df['Sequence Id']):
        df.at[i, 'Sequence Id'] = contigs[chrom]
    gbk = {contig: {} for contig in sorted(set(df['Sequence Id']))}
    for i, row in df.iterrows():
        if pd.isna(row['Gene']):
            gene = 'Unknown'
        else:
            gene = row['Gene']
        if pd.isna(row['Product']):
            product = 'hypothetical protein'
        else:
            product = row['Product']
        gbk[row['Sequence Id']][(row['Start'], row['Stop'])] = (gene, product)

    out = pd.DataFrame(columns=['cluster', 'chrom', 'pos', 'gene', 'product'])
    cluster = pd.read_csv('clusters_'+abb+'.csv')
    for j, (i, c, p) in enumerate(zip(cluster['id'], cluster['chrom'], cluster['pos'])):
        a = annotate_pos(gbk, c, p)
        if a:
            pass
        else:
            a = ['Not annotated', 'Not annotated']
        out.loc[j] = [i, c, p, a[0], a[1]]
    out = out.drop_duplicates(subset=['cluster', 'product'])
    print(out.to_string(index=False))


def t_test():
    """T test for testing variant richness significance."""
    df = get_variants(join('..', 'variants', 'variants_comp_mapping.csv'), 0)
    at = df[df['strain'] == strains[s.abbreviations['at']]]
    for j in ['T11', 'T22', 'T33', 'T44']:
        g1 = at[(at['timepoint'] == j) & (
            at['treatment'] == 3)]['hill'].to_numpy()
        g2 = at[(at['timepoint'] == j) & (
            at['treatment'] == 4)]['hill'].to_numpy()
        t, p = ttest_ind(g1, g2)
        if p < 0.05:
            print('At', j)
            print('T-Statistic', t, 'P-Value', p)

    ct = df[df['strain'] == strains[s.abbreviations['ct']]]
    for j in ['T11', 'T22', 'T33', 'T44']:
        g1 = ct[(ct['timepoint'] == j) & (
            ct['treatment'] == 2)]['hill'].to_numpy()
        g2 = ct[(ct['timepoint'] == j) & (
            ct['treatment'] == 3)]['hill'].to_numpy()
        g3 = ct[(ct['timepoint'] == j) & (
            ct['treatment'] == 4)]['hill'].to_numpy()
        t, p = ttest_ind(g1, g2)
        if p < 0.05:
            print('Ct treatment 2 and 3', j)
            print('T-Statistic', t, 'P-Value', p)
        t, p = ttest_ind(g1, g3)
        if p < 0.05:
            print('Ct treatment 2 and 4', j)
            print('T-Statistic', t, 'P-Value', p)
        t, p = ttest_ind(g2, g3)
        if p < 0.05:
            print('Ct treatment 3 and 4', j)
            print('T-Statistic', t, 'P-Value', p)


def plotter():
    variants = join('..','variants','variants_comp_mapping.csv')
    snps = join('..','variants','snps_comp_mapping.csv')
    fig = diversity_illumina(snps, 'Fixed variant richness','')
    fig = diversity_illumina(variants, 'Variant richness','')
    fig = snp_distribution(variants,'ct','T44',add_clusters=True)
    fig = snp_distribution(variants,'at','T44',add_clusters=True)
    fig = coverage()
    fig = trajectories('at','')
    fig = trajectories('ct','')
    fig = ct_box(variants,'Variants','')
    fig = ct_box(snps,'Variants','')

plotter()