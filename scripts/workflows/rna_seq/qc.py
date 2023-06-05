import plotly.express as px
import pysam
from samples import Samples
import pandas as pd
from os.path import split, join
from Bio import SeqIO
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



def annnotate(genes):
    products = []
    genes = ['NHFMJM_' + g for g in genes]
    f = 'ct.assembly.contigs.polypolish.transcript.gtf'
    columns = ['contig', 'source', 'type', 'start',
            'stop', 'a', 'strand', 'b', 'product']

    df = pd.read_csv(f, names=columns, sep='\t')

    f = 'assembly.contigs.polypolish.gbff'
    features = {}
    gbk = [contig for contig in SeqIO.parse(f, 'genbank')]
    for contig in gbk:
        for feature in contig.features:
            try:
                features[feature.qualifiers['locus_tag'][0]
                        ] = feature.qualifiers['product']
            except KeyError:
                pass

    for p in df['product']:
        for g in genes:
            if g in p:
                print(g,features[g][0])
                products.append(features[g][0])
    return products

Ct4 = annnotate(gs_ext)
Ct2 = annnotate(ct2)
