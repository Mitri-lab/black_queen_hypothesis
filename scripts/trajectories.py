import vcfpy
from samples import Samples
from os.path import join, exists
import pandas as pd
import plotly.express as px
from Bio import SeqIO

s = Samples()


class SNPs:
    def __init__(self, filter=20):
        self.snps = {strain: [] for strain in s.strains}
        for strain, samples in s.strains.items():
            for sample in samples:
                if sample['platform'] == 'illumina':
                    snps = []
                    f = join(sample['dir_name'], 'var.vcf')
                    if exists(f):
                        reader = vcfpy.Reader.from_path(f)
                        for record in reader:
                            snp = {}
                            snp['chrom'] = record.CHROM
                            snp['pos'] = record.POS
                            snp['qual'] = record.QUAL
                            snp['depth'] = record.INFO['DP']
                            snp['freq_sum'] = sum(
                                record.INFO['AO'])/record.INFO['DP']
                            snp['alt_depth_sum'] = sum(record.INFO['AO'])
                            snp['alt_depth'] = record.INFO['AO']
                            snp['freq'] = [freq / record.INFO['DP']
                                           for freq in record.INFO['AO']]
                            snp['alt'] = record.ALT
                            if filter is None:
                                snps.append(snp)
                            else:
                                if snp['qual'] >= filter:
                                    snps.append(snp)
                    else:
                        if strain == s.abbreviations['at']:
                            print(f)
                    sample['snps'] = snps
                    self.snps[strain].append(sample)
    
    def join_genbank(self):
        genbanks = {}
        for strain in s.strains:
            f = s.references[strain].replace('.fasta','_stripped.gbk')
            contigs = [contig for contig in SeqIO.parse(
            f, 'genbank')]
            genbank = {contig.id: {} for contig in contigs}
            for contig in contigs:
                for feature in contig.features:
                    # Some features don't have all desired keys
                    try:
                        start = feature.location.start
                        end = feature.location.end
                        product = feature.qualifiers['product']
                        genbank[contig.id][(start, end)] = product[0]
                    except KeyError:
                        pass
            genbanks[strain] = genbank
        for strain, samples in self.snps.items():
            for sample in samples:
                pass
        return genbanks


    def plot_phred_score(self):
        df = pd.DataFrame(list(), columns=['cosm', 'qual', ])
        """for i in df.index:
            df.at[i,'qual'] = []"""
        i = 0
        for sample in self.snps[s.abbreviations['at']]:
            for snp in sample['snps']:
                df.at[i, 'qual'] = snp['qual']
                df.at[i, 'cosm'] = str(
                    sample['treatment']) + '.' + str(sample['cosm'])
                i += 1
        fig = px.histogram(df, x='qual', facet_col='cosm', log_y=True,
                           facet_col_wrap=4, title='SNP phred scores')
        fig.write_image('hist_20.png', scale=2)

    def plot_trajectories(self):
        df = pd.DataFrame(columns=['micro_treat', 'timepoint', 'freq','qual','color'])
        i = 0
        for sample in self.snps[s.abbreviations['at']]:
            for snp in sample['snps']:
                key = '.'.join([str(item) for item in [snp['chrom'],
                            snp['pos'], sample['treatment'], sample['cosm']]])
                df.at[i, 'freq'] = snp['freq_sum']
                df.at[i, 'micro_treat'] = str(
                    sample['treatment']) + '.' + str(sample['cosm'])
                df.at[i, 'timepoint'] = sample['timepoint']
                df.at[i,'qual'] = snp['qual']
                df.at[i,'color'] = key
                i += 1
        mask = []
        for i in df['color']:
            if 'AGTU001.0001.c01.2351' in i:
                mask.append(True)
            else:
                mask.append(False)
        #df = df[mask]
        fig = px.line(df, x='timepoint', y='freq',color='color', facet_col='micro_treat',facet_col_wrap=3)
        fig.update_xaxes(categoryorder='array', categoryarray= ['T11','T22','T33','T44'])
        return df

    

snps = SNPs(filter=None)
# snps.plot_phred_score()
"""freqs = snps.plot_trajectories()
cosms = set(freqs['micro_treat'])
for cosm in cosms:
    tmp = freqs[freqs['micro_treat'] == cosm]
    print(cosm,len(tmp))"""
