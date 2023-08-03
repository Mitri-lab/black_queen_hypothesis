import pandas as pd
from os.path import join
from samples import Samples

s = Samples()

def annotate_clusters(abb):
    """Prints annotations for identified SNPs."""
    def annotate_pos(gbk, c, p):
        for (start, end), (gene, product) in gbk[c].items():
            if p in range(start, end):
                return [gene, product]
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

    out = pd.DataFrame(columns=['chrom', 'pos', 'gene', 'product','treatment','timepoint','cosm'])
    cluster = pd.read_csv(join('..', 'variants', 'snps_freebayes_comp_mapping.csv'))
    cluster = cluster[cluster['strain'] == s.abbreviations[abb]]
    for i,row in cluster.iterrows():
        a = annotate_pos(gbk, row['chrom'],row['pos'])
        if a:
            pass
        else:
            a = ['Not annotated', 'Not annotated']
        out.loc[i] = [row['chrom'],row['pos'],a[0],a[1],row['treatment'],row['timepoint'],row['cosm']]
    out.to_csv(join('..','annotations',abb+'_snps_annotations.csv'),index=False)

f = join('..','annotations','at_snps_annotations.csv')
df = pd.read_csv(f)
df = df[df['timepoint'] == 'T44']
products = set(df['product'])
out = pd.DataFrame(columns=['Gene','Product','Treatments'])
for p in products:
    tmp = df[df['product'] == p]
    if len(tmp) > 1:
        tmp.index = range(len(tmp))
        ts = ', '.join([str(e) for e in (list(set(tmp['treatment'])))])
        out.loc[len(out)] = [tmp.loc[0]['gene'],p,ts]