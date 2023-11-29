import pandas as pd
from os.path import join
from samples import Samples
from Bio import SeqIO
import plotly.express as px
import numpy as np
s = Samples()

w = 400
h = 250


def annotate_variants(abb):
    """Prints annotations for identified SNPs."""
    def annotate_pos(gbk, c, p):
        for (start, end), (gene, product) in gbk[c].items():
            if p in range(start, end):
                return [(start, end), gene, product]
        return False

    f = join('..', 'annotations', abb + '_renamed.tsv')
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

    in_files = [join('..', 'variants', 'variants_comp_mapping.csv'), join(
        '..', 'variants', 'snps_freebayes_comp_mapping.csv'),
        join('..', 'variants', 'snps_pacbio_renamed.csv')]
    out_files = [join('..', 'annotations', abb + '_variants_annotations.csv'),
                 join('..', 'annotations', abb + '_snps_annotations.csv'),
                 join('..', 'annotations', abb + '_pacbio_annotations.csv')]
    for j, in_file in enumerate(in_files):
        df = pd.read_csv(in_file)
        df = df[df['strain'] == s.abbreviations[abb]]
        df.insert(len(df.columns), 'gene_start', None)
        df.insert(len(df.columns), 'gene_end', None)
        df.insert(len(df.columns), 'gene', None)
        df.insert(len(df.columns), 'product', None)
        df.insert(len(df.columns), 'sequence', None)
        for i, row in df.iterrows():
            a = annotate_pos(gbk, row['chrom'], row['pos'])
            if a:
                pass
            else:
                a = [('Not annotated', 'Not annotated'),
                     'Not annotated', 'Not annotated']
            df.at[i, 'gene_start'], df.at[i, 'gene_end'] = a[0][0], a[0][1]
            df.at[i, 'gene'], df.at[i, 'product'] = a[1], a[2]
        contigs = {contig.name: str(contig.seq) for contig in SeqIO.parse(
            join('..', 'annotations', abb+'.fasta'), 'fasta')}
        for i, row in df.iterrows():
            if row['gene'] != 'Not annotated':
                seq = contigs[row['chrom']][int(
                    row['gene_start']):int(row['gene_end'])]
            else:
                seq = 'Not annotated'
            df.at[i, 'sequence'] = seq
        df.to_csv(out_files[j], index=False)


def filter():
    fs = [join('..', 'annotations', 'at_snps_annotations.csv'),
          join('..', 'annotations', 'ct_snps_annotations.csv')]
    outs = []
    for f in fs:
        print(f)
        df = pd.read_csv(f)
        df = df[df['timepoint'] == 'T44']
        products = set(df['product'])
        out = pd.DataFrame(
            columns=['Gene', 'Start', 'End', 'freq', 'Count', 'Product', 'Eff', 'Treatments', 'Sequence'])
        for p in products:
            tmp = df[df['product'] == p]
            if len(tmp) > 0:
                tmp.index = range(len(tmp))
                ts = ', '.join([str(e) for e in (list(set(tmp['treatment'])))])
                counts = len(tmp)
                out.loc[len(out)] = [tmp.loc[0]['gene'], tmp.loc[0]['start'],
                                     tmp.loc[0]['end'], tmp.loc[0]['freq'], counts, p, tmp.loc[0]['eff'], ts, tmp.loc[0]['sequence']]
        outs.append(out)
    return outs


def tex_table():
    outs = filter()
    catpions = ['Fixed mutations for At', 'Fixed mutations for Ct']
    for i, o in enumerate(outs):
        o = o[['Gene', 'Product', 'Treatments', 'Eff']]
        o.columns = ['Gene', 'Product', 'Treatment', 'Effect']
        print(o.to_latex(index=False, caption=catpions[i]))


def font_size(fig):
    """Style function for figures setting fot size and true black color."""
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


def products(abb):
    df = pd.read_csv(join('..', 'annotations', abb+'_variants_annotations.csv'))
    filter = (df['timepoint'] == 'T44') & (
        df['freq'] >= 0.1) & (df['product'] != 'Not annotated')
    df = df[filter]
    products = sorted(
        list(set(zip(df['product'], df['gene']))), key=lambda x: x[0].lower())
    mult_freqs = []
    matrix = []
    sample_labels = []
    for sample in s.strains[s.abbreviations[abb]]:
        if (sample['timepoint'] == 'T44') & (sample['platform'] == 'illumina'):
            sample_labels.append(sample['name'])
            row = []
            tmp = df[df['name'] == sample['name']]
            for product, gene in products:
                filter = (tmp['product'] == product) & (tmp['gene'] == gene)
                i = tmp[filter]
                freqs = i['freq'].to_list()
                if len(freqs) > 0:
                    row.append(np.average(freqs))
                    if len(freqs) > 1:
                        mult_freqs.append((product, gene))
                else:
                    row.append(0)
            matrix.append(row)

    dm = pd.DataFrame(matrix)
    dm.index = [s.labels[key] for key in sample_labels]
    dm.columns = [' <i>' + gene + '</i> ' +
                  product for product, gene in products]

    colors = ['#7570B3', '#1B9E77', '#E6AB02']
    # colors = ['#7570B3', '#1B9E77', '#D95F02', '#E6AB02']

    for c in dm.columns:
        if len(dm[dm[c] != 0]) <= 1:
            dm = dm.drop(c, axis=1)
    custom_colorscale = ['white','black']
    fig = px.imshow(dm, color_continuous_scale=custom_colorscale)
    fig.update_traces(colorbar=dict(thickness=1))

    fig = font_size(fig)
    fig.update_layout(
        xaxis=dict(tickangle=-45),  # Rotate x-axis labels by -45 degrees
    )
    fig.write_image(join('..', 'plots', 'plots', abb + '_products.pdf'))
    fig.show()


k = []


def annotate_pacbio():
    contigs = {'tig00000001_polypolish': '_0',
               'tig00000002_polypolish': '_1',
               'tig00000003_polypolish': '_2',
               'tig00000004_polypolish': '_3',
               'tig00000005_polypolish': '_4'}

    mapping = {s.abbreviations['at']: 'at',
               s.abbreviations['ct']: 'ct',
               s.abbreviations['oa']: 'oa',
               s.abbreviations['ms']: 'ms'}

    df = pd.read_csv(join('..', 'variants', 'snps_pacbio.csv'))
    for i, (c, n) in enumerate(zip(df['chrom'], df['strain'])):
        df.at[i, 'chrom'] = mapping[n] + contigs[c]

    df.to_csv(join('..', 'variants', 'snps_pacbio_renamed.csv'), index=False)

    annotate_variants('ct')

    df = pd.read_csv(join('..', 'annotations', 'ct_pacbio_annotations.csv'))
    ct1 = df[df['name'] == 'Ct42.1']
    ct2 = df[df['name'] == 'Ct42.2']
