from samples import Samples
from os.path import join, exists
import vcfpy
import pandas as pd
import plotly.express as px

s = Samples()

strains = {s.abbreviations['at']: '<i>A. tumefaciens</i>',
           s.abbreviations['ct']: '<i>C. testosteroni</i>',
           s.abbreviations['oa']: '<i>O. anthropi</i>',
           s.abbreviations['ms']: '<i>M. saperdae</i>'}


def parse_variants():
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
                            snp['freq'] = record.INFO['AB'][i]
                            snp['alt'] = r
                            key = '.'.join([snp['chrom'], str(snp['pos']),
                                            str(sample['treatment']), str(sample['cosm']), str(snp['alt'])])
                            snp['linegroup'] = key
                            if snp['qual'] >= 20:
                                df = (pd.concat(
                                    [pd.DataFrame(snp, index=[0]), pd.DataFrame(sample, index=[0])], axis=1))
                                dfs.append(df)
    out = pd.concat(dfs)
    out = out[out['freq'] != 0]
    out['treatment'] = out['treatment'].astype(str)
    out['cosm'] = out['cosm'].astype(str)
    out.to_csv('variants.csv', index=False)


def depth():
    df = pd.DataFrame(columns=['sample', 'timepoint',
                               'species', 'treatment', 'depth', 'cosm'])
    i = 0
    for strain, samples in s.strains.items():
        for sample in samples:
            if sample['platform'] == 'illumina':
                f = join(sample['dir_name'], 'coverage.txt')
                with open(f, 'r') as handle:
                    d = float(handle.read())
                r = [sample['name'], sample['timepoint'], strains[strain],
                     sample['treatment'], d, sample['cosm']]
                df.loc[i] = r
                i += 1
    df.to_csv(
        '/users/eulrich/work/genome_size/data/bulk/depth.csv', index=False)


def diversity(q):
    snps = pd.read_csv('variants.csv')
    titles = ['T11', 'T22', 'T33', 'T44']
    hill = pd.DataFrame(
        columns=['strain', 'treatment', 'hill', 'timepoint', 'cosm'])
    i = 0
    for strain in s.strains.keys():
        tmp = snps[snps['strain'] == strain]
        for sample in set(tmp['name']):
            sub = tmp[tmp['name'] == sample]
            j = 0
            if q == 0:
                for p, t, f, c in zip(sub['timepoint'], sub['treatment'], sub['freq'], sub['cosm']):
                    j += 1
                try:
                    hill.loc[i] = [strains[strain], t, j, p, c]
                    i += 1
                except ZeroDivisionError:
                    pass

            elif q == 2:
                for p, t, f, c in zip(sub['timepoint'], sub['treatment'], sub['freq'], sub['cosm']):
                    j += f**2
                try:
                    hill.loc[i] = [strains[strain], t, 1/j, p, c]
                    i += 1
                except ZeroDivisionError:
                    pass
    hill = hill.sort_values(by='treatment', ascending=True)
    hover_data = ['treatment', 'timepoint', 'strain', 'hill', 'cosm']
    n_colors = len(set(hill['treatment']))
    colors = px.colors.sample_colorscale(
        "Agsunset", [n/(n_colors - 1) for n in range(n_colors)])
    fig = px.strip(hill, x='strain', y='hill', color='treatment', color_discrete_sequence=colors,
                   stripmode="overlay", facet_col='timepoint',
                   log_y=False, category_orders={'timepoint': ['T11', 'T22', 'T33', 'T44']},
                   hover_data=hover_data)
    for d in fig['data']:
        d['marker']['size'] = 8
    for i, t in enumerate(fig['layout']['annotations']):
        t['text'] = titles[i]

    fig.for_each_yaxis(lambda y: y.update(title=''))
    fig['layout']['yaxis']['title']['text'] = 'N variants'
    fig.update_layout(font={'size': 14})
    fig.update_xaxes(title=None)
    fig.show()


def heterogenity(abb, timepoint, subset=False, add_clusters=False):
    snps = pd.read_csv('variants.csv', dtype={'cosm': str, 'treatment': str})
    strain = snps[snps['strain'] == s.abbreviations[abb]]
    strain = strain[strain['timepoint'] == timepoint]
    contigs = {j: 'Contig ' +
               str(i) for i, j in enumerate(sorted(list(set(strain['chrom']))))}
    strain.index = range(len(strain))
    strain.insert(0, 'ID', strain.index)
    for i, c in enumerate(strain['chrom']):
        strain.at[i, 'chrom'] = contigs[c]
    if subset:
        strain = strain[strain['treatment'] == subset]
    hover_data = ['pos', 'depth', 'qual', 'cosm', 'ID']

    if subset:
        n_colors = len(set(strain['cosm']))
        colors = px.colors.sample_colorscale(
            "Agsunset", [n/(n_colors - 1) for n in range(n_colors)])
        fig = px.scatter(strain, x='pos', y='freq', color='cosm', color_discrete_sequence=colors,
                         facet_col='chrom', facet_col_wrap=1, category_orders={'chrom': contigs.values()},
                         facet_row_spacing=0.12, hover_data=hover_data)
    else:
        n_colors = len(set(strain['cosm']))
        colors = px.colors.sample_colorscale(
            "Agsunset", [n/(n_colors - 1) for n in range(n_colors)])
        fig = px.scatter(strain, x='pos', y='freq', color='treatment', color_discrete_sequence=colors,
                         facet_col='chrom', facet_col_wrap=1, category_orders={'chrom': contigs.values()},
                         facet_row_spacing=0.12, hover_data=hover_data)
    fig.update_xaxes(matches=None, showticklabels=True)
    fig.for_each_yaxis(lambda y: y.update(title=''))
    fig.add_annotation(x=-0.05, y=0.5,
                       text="Variant frequency", textangle=-90,
                       xref="paper", yref="paper")
    fig.for_each_xaxis(lambda x: x.update(title=''))
    fig.add_annotation(x=0.5, y=-0.1,
                       text="Position", textangle=0,
                       xref="paper", yref="paper", showarrow=False)
    for i, t in enumerate(fig['layout']['annotations']):
        t['text'] = t['text'].replace('chrom=', '')
    if subset:
        for d in fig['data']:
            d['marker']['size'] = 6
    else:
        for d in fig['data']:
            d['marker']['size'] = 6
    fig.update_layout(font={'size': 14})
    if add_clusters:
        annot = pd.read_csv('clusters_'+abb+'.csv')
        for i, chrom in enumerate(sorted(set(annot['chrom']))):
            annot_sub = annot[annot['chrom'] == chrom]
            for id in set(annot_sub['id']):
                if (abb == 'at') & (chrom == 'tig00000002_polypolish'):
                    i = 3
                asid = annot_sub[annot_sub['id'] == id]
                asid.index = range(len(asid))
                x0 = asid.loc[0].pos
                x1 = asid.loc[len(asid)-1].pos
                fig.add_vrect(x0=x0, x1=x1, row=i, annotation_text=id,
                              opacity=0.25)
    fig.show()


def coverage():
    df = pd.read_csv('depth.csv')
    df = df.sort_values(by='treatment', ascending=True)
    hover_data = ['treatment', 'timepoint', 'species', 'depth', 'cosm']
    n_colors = len(set(df['treatment']))
    colors = px.colors.sample_colorscale(
        "Agsunset", [n/(n_colors - 1) for n in range(n_colors)])
    fig = px.strip(df, x='species', y='depth', color='treatment', facet_col='timepoint',
                   category_orders={'timepoint': [
                       'T11', 'T22', 'T33', 'T44']}, color_discrete_sequence=colors,
                   log_y=True, hover_data=hover_data)
    for d in fig['data']:
        d['marker']['size'] = 8
    for i, t in enumerate(fig['layout']['annotations']):
        t['text'] = titles[i]

    fig.for_each_yaxis(lambda y: y.update(title=''))
    fig['layout']['yaxis']['title']['text'] = 'Mean coverage'
    fig.update_layout(font={'size': 14})
    fig.update_xaxes(title=None)
    fig.show()


def trajectories(species):
    df = pd.read_csv('variants.csv', dtype={'cosm': str, 'treatment': str})
    df = df[df['strain'] == s.abbreviations[species]]
    n_colors = len(set(df['cosm']))
    colors = px.colors.sample_colorscale(
        "Agsunset", [n/(n_colors - 1) for n in range(n_colors)])
    hover_data = ['treatment', 'timepoint', 'cosm', 'depth']
    fig = px.line(df, x='timepoint', y='freq', line_group='linegroup', color_discrete_sequence=colors,
                  facet_col='treatment', facet_col_wrap=4, color='cosm', hover_data=hover_data,
                  category_orders={'timepoint': ['T11', 'T22', 'T33', 'T44']}, markers=True)
    for d in fig['data']:
        d['marker']['size'] = 3
        d['line']['width'] = 1
    fig.for_each_yaxis(lambda y: y.update(title=''))
    fig['layout']['yaxis']['title']['text'] = 'Variant frequency'
    fig.update_xaxes(title=None)
    fig.update_layout(font={'size': 14})
    titles = ['Treatment ' + str(i)
              for i in sorted(list(set(df['treatment'])))]
    for i, t in enumerate(fig['layout']['annotations']):
        t['text'] = titles[i]
    fig.add_annotation(x=0.5, y=-0.1,
                       text="Timepoint", textangle=0,
                       xref="paper", yref="paper", showarrow=False)
    fig.show()


def cluster_trajectories(abb):
    df = pd.read_csv('variants.csv', dtype={'cosm': str, 'treatment': str})
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
                  facet_col='treatment', facet_col_wrap=4, color='cluster', color_discrete_sequence=colors,
                  category_orders={'timepoint': ['T11', 'T22', 'T33', 'T44'], 'treatment': sorted(set(out['treatment']))}, markers=True)
    titles = ['Treatment ' + str(i)
              for i in sorted(list(set(df['treatment'])))]
    for i, t in enumerate(fig['layout']['annotations']):
        t['text'] = titles[i]
    fig.for_each_yaxis(lambda y: y.update(title=''))
    fig['layout']['yaxis']['title']['text'] = 'Variant frequency'
    fig.update_xaxes(title=None)
    fig.update_layout(font={'size': 14})
    for d in fig['data']:
        d['marker']['size'] = 3
        d['line']['width'] = 1
    fig.add_annotation(x=0.5, y=-0.1,
                       text="Timepoint", textangle=0,
                       xref="paper", yref="paper", showarrow=False)
    fig.show()


def annotate_clusters(abb):
    def annotate_pos(gbk, c, p):
        for (start, end), (gene, product) in gbk[c].items():
            if p in range(start, end):
                return [gene, product]
        return False

    f = abb + '.tsv'
    df = pd.read_csv(f, sep='\t')
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
