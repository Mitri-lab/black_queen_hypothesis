from experiment import Experiment
from samples import Samples
from os.path import join
from os.path import exists
import pandas as pd

e = Experiment()
s = Samples()

class Mutations():
    def __init__(self,dropna=False):
        self.strains = {key:None for key in s.strains.keys()}
        for strain,samples in s.strains.items():
            dfs = []
            for sample in samples:
                f = join(sample['dir_name'],'snippy','snps.tab')
                if exists(f):
                    df = pd.read_csv(f,sep='\t')
                    if dropna:
                        df = df.dropna()
                    df.insert(0,'sample',sample['name'])
                    dfs.append(df)
            self.strains[strain] = pd.concat(dfs)

    def get_gene_series(self):
        """Returns a dataframe per treatment and strain
        with timepoints as x axis and genes as y axis,
        values are observations in n cosms."""
        self.gene_series = {key:None for key in e.strain_series.keys()}
        #First we find all affected genes in per series
        #and create out DataFrame
        for (strain,treatment),samples in e.strain_series.items():
            dfs = []
            for sample in samples:
                f = join(sample['dir_name'],'snippy','snps.tab')
                if exists(f):
                    dfs.append(pd.read_csv(f,sep='\t'))
            genes = set(pd.concat(dfs).dropna()['GENE'])
            out = pd.DataFrame(0,columns=e.timepoints,index=genes)
            #Filling created dataframe
            for sample in samples:
                f = join(sample['dir_name'],'snippy','snps.tab')
                if exists(f):
                    genes = pd.read_csv(f,sep='\t').dropna()['GENE']
                    for gene in genes:
                        out.at[gene,sample['timepoint']] += 1
            self.gene_series[(strain,treatment)] = out