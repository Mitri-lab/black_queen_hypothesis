import plotly.express as px
import pysam
from samples import Samples
import pandas as pd
from os.path import split, join
from Bio import SeqIO
import math
s = Samples()


def dump_depth():
    df = pd.DataFrame(columns=['sample', 'depth'])
    j = 0
    for strain, samples in s.strains.items():
        for sample in samples:
            if sample['platform'] == 'RNA':
                print(sample['name'])
                i = 0
                a = pysam.AlignmentFile(sample['dir_name'])
                for read in a:
                    if not (read.is_unmapped):
                        i += 1
                df.loc[j] = [sample['name'], str(i)]
                j += 1
    df.to_csv(join(split(sample['dir_name'])[0], 'coverage.csv'), index=False)


def plot_depth():
    df = pd.read_csv('depth.csv')
    df = df.sort_values(by='sample', ascending=True)
    fig = px.scatter(df, x='sample', y='depth', log_y=True)
    fig.show()


gs = ['21645', '23650', '12335', '08300', '08305', '14200',
      '00170', '14195', '06560', '18990', '19650', '19805',
      '02180', '06555', '09295', '01885', '19810', '14185',
      '19810', '14185', '14190', '02175']

gs_ext = ['19650','06555','14185','23650','14190','08305','12335']

ct2_ext = ['03685','20235','06455','11590','24860','02180','25795','04460','06555','06560','01885','18275','02700','04915','21645','20240','10260','13355','16850','14185']
ct2 = ['02180','20240','06555','20235','18275','06560','04460','10260']


def annotate(names):
    products = []
    genes = ['NHFMJM_' + g for g in names]
    f = 'ct.assembly.contigs.polypolish.transcript.gtf'
    columns = ['contig', 'source', 'type', 'start',
            'stop', 'a', 'strand', 'b', 'product']

    df = pd.read_csv(f, names=columns, sep='\t')

    f = 'assembly.contigs.polypolish.gbff'
    features = {}
    gbk = [contig for contig in SeqIO.parse(f, 'genbank')]
    for contig in gbk:
        for feature in contig.features:
            entry = [None,None]
            try:
                entry[0] = feature.qualifiers['product'][0]
            except KeyError:
                entry[0] = 'Unknown'
            try:
                entry[1] = feature.qualifiers['gene'][0]
            except KeyError:
                entry[1] = 'Unknown'
            try:
                features[feature.qualifiers['locus_tag'][0]] = entry
            except KeyError:
                pass
    annot = pd.read_csv('deg_CtA_Ct2.tsv', sep='\t')
    annot.index = annot['genes']
    for p in df['product']:
        for g in genes:
            if g in p:
                fc = annot.loc[g]['logFC']
                pv = -math.log10(annot.loc[g]['FDR'])
                entry = (g,features[g][1],features[g][0],fc,pv)
                products.append(entry)
    out = pd.DataFrame(columns=['Gene ID','Gene','Product','Log2 fold change','-Log10P'])
    for i in products:
        out.loc[len(out)] = list(i)
    return out

df = annotate(ct2)
caption = 'Differentially expressed genes in mono-evolved \textit{C. testosteroni}'
print(df.to_latex(index=False,caption=caption))
