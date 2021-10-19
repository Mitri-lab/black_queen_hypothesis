from samples import Samples
import pandas as pd
from os.path import join
from os.path import exists

s = Samples()

class Mutations():
    def get_snps(self):
        self.strains = {strain:None for strain in s.strains.keys()}
        for strain in s.strains:
            snps = Snps(strain)
            snps.get_avg_snps()
            self.strains[strain] = snps

class Snps():
    """This class returns SNPs stats per strain."""
    def __init__(self,strain):
        self.strain = strain
        #Grabs all treatments per strain
        self.treatments = set([sample['treatment'] for sample in s.strains[self.strain]])

        #Stors all snippy dfs of a strain with the treatment as a key
        self.snippy_dfs = {treatment:[] for treatment in self.treatments}
        for sample in s.strains[self.strain]:
            #We only look at Illumina and samples with a snippy file
            #(because of no coverage some samples have no snippy file)
            if (sample['platform'] == 'illumina') & exists(join(sample['dir'],'snippy','snps.tab')):
                df = pd.read_csv(join(sample['dir'],'snippy','snps.tab'),sep='\t')
                self.snippy_dfs[sample['treatment']].append(df)

        #Stors n_samples per treatment for average calculations
        self.samples_per_treatment = {treatment:None for treatment in self.treatments}
        for treatment in self.treatments:
            n_samples = len([sample for sample in s.strains[self.strain] \
            if (sample['treatment'] == treatment) & (sample['platform'] == 'illumina')])
            self.samples_per_treatment[treatment] = n_samples


    def get_avg_snps(self):
        """Calculates average SNPS per treatment"""
        self.avg_snps = dict()
        n_snps = dict()
    
        #Setting values of helper dictionaries
        for treatment in self.treatments:
            n_snps[treatment] = []
            self.avg_snps[treatment] = None
        
        #Calculating average snps per treatment
        for treatment in self.snippy_dfs:
            for df in self.snippy_dfs[treatment]:
                n_snps[treatment].append(len(df))
            n_samples = self.samples_per_treatment[treatment]
            self.avg_snps[treatment] = sum(n_snps[treatment])/n_samples

    def get_avg_effects(self,strain):
        pass

m = Mutations()
m.get_snps()