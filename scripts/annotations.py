import pandas as pd
from os.path import join
from samples import Samples
from Bio import SeqIO
s = Samples()


def annotate_clusters(abb):
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

    out = pd.DataFrame(columns=['chrom', 'start', 'end',
                       'gene', 'product', 'treatment', 'timepoint', 'cosm'])
    cluster = pd.read_csv(
        join('..', 'variants', 'snps_freebayes_comp_mapping.csv'))
    cluster = cluster[cluster['strain'] == s.abbreviations[abb]]
    for i, row in cluster.iterrows():
        a = annotate_pos(gbk, row['chrom'], row['pos'])
        if a:
            pass
        else:
            a = [('Not annotated', 'Not annotated'),
                 'Not annotated', 'Not annotated']
        out.loc[i] = [row['chrom'], a[0][0], a[0][1], a[1], a[2],
                      row['treatment'], row['timepoint'], row['cosm']]
    out.insert(len(out.columns), 'sequence', None)
    contigs = {contig.name: str(contig.seq) for contig in SeqIO.parse(
        join('..', 'annotations', abb+'.fasta'), 'fasta')}
    for i, row in out.iterrows():
        if row['gene'] != 'Not annotated':
            seq = contigs[row['chrom']][int(row['start']):int(row['end'])]
        else:
            seq = 'Not annotated'
        out.at[i, 'sequence'] = seq
    out.to_csv(join('..', 'annotations', abb +
               '_snps_annotations.csv'), index=False)
    return out

def selection():
    fs = [join('..', 'annotations', 'at_snps_annotations.csv'),
        join('..', 'annotations', 'ct_snps_annotations.csv')]
    outs = []
    for f in fs:
        print(f)
        df = pd.read_csv(f)
        df = df[df['timepoint'] == 'T44']
        products = set(df['product'])
        out = pd.DataFrame(
            columns=['Gene', 'Start', 'End','Count', 'Product', 'Treatments', 'Sequence'])
        for p in products:
            tmp = df[df['product'] == p]
            if len(tmp) > 1:
                tmp.index = range(len(tmp))
                ts = ', '.join([str(e) for e in (list(set(tmp['treatment'])))])
                counts = len(tmp)
                out.loc[len(out)] = [tmp.loc[0]['gene'], tmp.loc[0]['start'],
                                    tmp.loc[0]['end'],counts, p, ts, tmp.loc[0]['sequence']]
        outs.append(out)
    return outs


