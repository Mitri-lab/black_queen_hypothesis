from dna_features_viewer import GraphicFeature, GraphicRecord,CircularGraphicRecord
import matplotlib.pyplot as plt
import pandas as pd
from samples import Samples
from Bio import SeqIO
from os.path import join

s = Samples()

plt.rcParams["font.size"] = 10




def plot_genbank(genbank_list, chromosome, start,end, out):
    """Plots genbank annotation of deleted sequence"""
    genbank = {contig.id: contig for contig in genbank_list}
    features = []
    for feature in genbank[chromosome].features:
        if feature.type == 'gene':
            f_start = feature.location.start
            f_end = feature.location.end
            if (f_start >= start) & (f_end <= end):
                try:
                    gf = GraphicFeature(start=f_start, end=f_end, strand=feature.location.strand,
                                        color="#000000", label=feature.qualifiers['gene'][0])
                    features.append(gf)
                except KeyError:
                    pass
    record = GraphicRecord(sequence_length=end, features=features)
    record = record.crop((start, end))
    f_name = '.'.join([chromosome, str(start), 'pdf'])
    fig,axes = record.plot(figure_width=10,figure_height=5)
    return fig,axes
    """record.plot_on_multiple_pages(join(out, f_name),
                                  nucl_per_line=(end-start),
                                  lines_per_page=10,
                                  plot_sequence=False
                                  )
"""
def plot_deletion():
    pass

genbank = [contig for contig in SeqIO.parse('/users/eulrich/work/genome_size/data/ancestors/ct/bakta/assembly.contigs.polypolish.gbff', 'genbank')]

fig,ax = plot_genbank(genbank,'tig00000001_polypolish',5638639,5638639+145276,join('..','plots','plots'))

fig.figure.savefig('tmp2.svg')
"""features = []
for c in refs[abb].values():
    features.append(GraphicFeature(start=c[0], end=c[1], color="#ffd700"))

df = pd.read_csv('../variants/depth_Q_0.concat.csv')
for sample in s.strains[s.abbreviations[abb]]:
    if sample['platform'] == 'pacbio':
        df_sample = df[(df['name'] == sample['name']) & (df['length'] > 10)]
        for c in refs[abb]:
            tmp = df_sample[df_sample['contig'] == c]
            prev = refs[abb][c][0]
            c_0 = refs[abb][c][0]
            for start,end in zip(tmp['start'],tmp['end']):
                features.append(GraphicFeature(start=prev, end=start+c_0, color="#ffd700"))
                prev = end+c_0
            features.append(GraphicFeature(start=start+c_0, end=refs[abb][c][1], color="#ffd700"))



record = CircularGraphicRecord(sequence_length=list(refs[abb].values())[-1][1], features=features)
record.plot(figure_width=5)"""