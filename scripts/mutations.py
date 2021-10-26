from experiment import Experiment
from samples import Samples
from os.path import join
from os.path import exists
import pandas as pd

e = Experiment()
s = Samples()

class Mutations():
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
                    #If gene is mutated multiple times we want to count it only once
                    gene_added = {gene:False for gene in set(genes)}
                    for gene in genes:
                        if not gene_added[gene]:
                            out.at[gene,sample['timepoint']] += 1
                            gene_added[gene] = True
            self.gene_series[(strain,treatment)] = out