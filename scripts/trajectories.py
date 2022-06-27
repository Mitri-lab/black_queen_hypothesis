import vcfpy
from samples import Samples
from os.path import join, exists
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from Bio import SeqIO
from colour import Color

s = Samples()


class SNPs:
    def __init__(self, filter=20, parse=False):
        if parse:
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
                                snp['af'] = record.INFO['AF']
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
        else:
            self.df = pd.read_csv(join(s.work, 'bulk', 'snps.csv'))

    def join_genbank(self):
        genbanks = {}
        for strain in s.strains:
            f = s.references[strain].replace('.fasta', '_stripped.gbk')
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
                for snp in sample['snps']:
                    snp['product'] = 'intergenic'
                    for (start, end), product in genbanks[strain][snp['chrom']].items():
                        if (snp['pos'] >= start) & (snp['pos'] <= end):
                            snp['product'] = product

    def get_data_frame(self, strain, filter_genes=False):
        self.df = pd.DataFrame(
            columns=['micro_treat', 'timepoint', 'freq', 'qual', 'color', 'product'])
        i = 0
        for sample in self.snps[strain]:
            for snp in sample['snps']:
                key = '.'.join([str(item) for item in [snp['chrom'],
                                                       snp['pos'], sample['treatment'], sample['cosm']]])
                self.df.at[i, 'freq'] = snp['freq_sum']
                self.df.at[i, 'micro_treat'] = str(
                    sample['treatment']) + '.' + str(sample['cosm'])
                self.df.at[i, 'timepoint'] = sample['timepoint']
                self.df.at[i, 'qual'] = snp['qual']
                self.df.at[i, 'color'] = key
                self.df.at[i, 'product'] = snp['product']
                self.df.at[i, 'treatment'] = sample['treatment']
                i += 1

        if filter_genes:
            self.df = self.df.drop_duplicates(
                ['product', 'timepoint', 'micro_treat'])

    def plot_phred_score(self, out):
        fig = px.histogram(self.df, x='qual', facet_col='micro_treat', log_y=True,
                           facet_col_wrap=5, title='SNP phred scores',
                           category_orders={'micro_treat': list(sorted(self.df['micro_treat']))})
        fig.write_html(out + '_phred_score.html')

    def plot_trajectories(self, out):
        mask = []
        for snp in self.df['color']:
            if len(set(self.df[self.df['color'] == snp]['timepoint'])) > 1:
                mask.append(True)
            else:
                mask.append(False)
        df = self.df[mask]
        df = df.sort_values('product')

        products = {product: None for product in set(df['product'])}
        products = dict(sorted(products.items()))
        color_range = (Color("purple"), Color("orange"))
        spectrum = list(color_range[0].range_to(
            color_range[1], len(products.keys())))

        for product, color in zip(products, spectrum):
            products[product] = color.get_hex()

        colors = {}
        for snp, product in zip(df['color'], df['product']):
            colors[snp] = products[product]

        labels = {snp: product for snp, product in zip(
            df['color'], df['product'])}

        fig = px.line(df, x='timepoint', y='freq', line_group='color', color='product',
                      facet_col='micro_treat', facet_col_wrap=5,
                      hover_data=['product'],
                      category_orders={'micro_treat': list(
                          sorted(self.df['micro_treat']))},
                      color_discrete_map=products)
        fig.update_xaxes(categoryorder='array', categoryarray=[
                         'T11', 'T22', 'T33', 'T44'])
        # fig.update_traces(showlegend=True)
        fig.write_html(out + '_trajectories.html')
        """names = set()
        fig.for_each_trace(
            lambda trace:
                trace.update(showlegend=False)
                if (trace.name in names) else names.add(trace.name))"""
        return fig

    def plot_frequencies(self, out):
        fig = px.histogram(self.df, x='freq', facet_col='micro_treat', log_y=True,
                           facet_col_wrap=4, title='SNP frequencies',
                           category_orders={'micro_treat': list(sorted(self.df['micro_treat']))})
        fig.write_html(out + '_frequencies.html')

    def plot_freq_histogram(self, out):
        fig = px.histogram(self.df, x='freq', facet_col='micro_treat', log_y=False,
                           facet_col_wrap=5, title='Freq dist',
                           category_orders={'micro_treat': list(sorted(self.df['micro_treat']))})
        fig.write_html(out + '_frequencies_hist.html')

    def plot_freq_violin(self, species, out):
        treatments = s.treatments[species]
        fig = make_subplots(rows=1, cols=len(treatments),
                            shared_yaxes=True, subplot_titles=treatments)
        for counter, treatment in enumerate(treatments):
            for timepoint in ['T11', 'T22', 'T33', 'T44']:
                df = self.df[self.df['timepoint'] == timepoint]
                df = df[df['treatment'] == treatment]
                fig.add_trace(go.Violin(y=df['freq'], points='all', spanmode='hard', marker_color='blue',
                                        pointpos=0, name=timepoint, legendgroup=timepoint, showlegend=False), row=1, col=counter+1)
        fig.write_html(out + '_frequencies_violin.html')

    def plot_snps_bar(self, species, out):
        colors = {1: 'red',
                  2: 'red',
                  3: 'purple',
                  4: 'orange'}
        treatments = s.treatments[species]
        fig = go.Figure()
        for counter, treatment in enumerate(treatments):
            for timepoint in ['T11', 'T22', 'T33', 'T44']:
                df = self.df[self.df['timepoint'] == timepoint]
                df = df[df['treatment'] == treatment]
                n_samples = len(set(df['micro_treat']))
                fig.add_trace(go.Scatter(x=[timepoint], y=[len(
                    df)/n_samples], name=treatment, legendgroup=treatment, showlegend=True,line_color=colors[treatment], mode='markers'))
        fig.update_xaxes(type='category')
        fig.write_html(out + '_snps_bar.html')

    def plot_products(self, out):
        fig = px.histogram(self.df, x='product', facet_col='micro_treat', log_y=True,
                           facet_col_wrap=5, title='SNP frequencies',
                           category_orders={'micro_treat': list(sorted(self.df['micro_treat']))})
        fig.update_xaxes(showticklabels=False)
        fig.update_xaxes(categoryorder='array',
                         categoryarray=list(set(self.df['product'])))
        fig.write_html(out + '_products.html')


snps = SNPs(filter=20, parse=True)
snps.join_genbank()
strain = s.abbreviations['at']
out = s.abbreviations[strain]
snps.get_data_frame(strain)
snps.plot_freq_violin(strain, out)
snps.plot_snps_bar(strain, out)
# snps.plot_freq_histogram(out)
# snps.plot_products(out)
#fig = snps.plot_trajectories(out)
# snps.plot_phred_score(out)
# snps.plot_frequencies(out)
# snps.plot_trajectories(out)
# snps.get_data_frame(strain,filter_genes=True)
#out = 'unique_genes_' + out
# snps.plot_frequencies(out)
# snps.plot_trajectories(out)
"""freqs = snps.plot_trajectories()
cosms = set(freqs['micro_treat'])
for cosm in cosms:
    tmp = freqs[freqs['micro_treat'] == cosm]
    print(cosm,len(tmp))"""
