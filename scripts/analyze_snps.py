from samples import Samples
import pandas as pd
from os.path import join
from os.path import exists

s = Samples()

class Mutations():
    """Summarizes stats of SNPs over all strains."""
    def __init__(self):
        self.strains = {strain:None for strain in s.strains.keys()}

    def get_snps(self):
        for strain in s.strains:
            snps = Snps(strain)
            snps.get_avg_snps()
            snps.get_avg_effects()
            snps.get_gene_counts_per_treatment()
            self.strains[strain] = snps

    def print_effects(self):
        """Prints effects per treatments"""
        for strain,effects in self.strains.items():
            print(strain)
            print(pd.DataFrame(effects.avg_effects),'\n')

    def print_gene_counts(self):
        """Prints normalized gene counts per treatment."""
        for strain in self.strains:
            print(strain,'\n',pd.DataFrame(self.strains[strain].gene_counts_per_treatment))

class Snps():
    """This class returns SNPs stats per strain."""
    def __init__(self,strain):
        self.strain = strain

        #Timepoints of studies
        self.timepoints = ['T11','T22','T33','T44']
        #Grabs all treatments per strain
        self.treatments = set([sample['treatment'] for sample in s.strains[self.strain]])
        self.microcosms = [1,2,3,4,5]

        #Stores all snippy dfs of a strain with the treatment as a key
        #Stores all snippy dsf of a strain with the timepoint as a key
        #Stores all snippy dsf of a strain with the microcosm as a key
        self.snps_per_treatment = {treatment:[] for treatment in self.treatments}
        self.snps_per_timepoint = {timepoint:[] for timepoint in self.timepoints}
        self.snps_per_microcosm = {microcosm:[] for microcosm in self.microcosms}
        for sample in s.strains[self.strain]:
            #We only look at Illumina and samples with a snippy file
            #(because of no coverage some samples have no snippy file)
            if (sample['platform'] == 'illumina') & exists(join(sample['dir'],'snippy','snps.tab')):
                df = pd.read_csv(join(sample['dir'],'snippy','snps.tab'),sep='\t')
                self.snps_per_treatment[sample['treatment']].append(df)
                timepoint = sample['name'][0:3]
                self.snps_per_timepoint[timepoint].append(df)
                self.snps_per_microcosm[sample['microcosm']].append(df)

    def get_avg_snps(self):
        """Calculates average SNPS per treatment"""
        self.avg_snps = dict()
        n_snps = dict()
    
        #Setting values of helper dictionaries
        for treatment in self.treatments:
            n_snps[treatment] = []
            self.avg_snps[treatment] = None
        
        #Calculating average snps per treatment
        for treatment in self.snps_per_treatment:
            for df in self.snps_per_treatment[treatment]:
                n_snps[treatment].append(len(df))
            n_samples = len(self.snps_per_treatment[treatment])
            self.avg_snps[treatment] = sum(n_snps[treatment])/n_samples

    def get_avg_effects(self):
        """Counts effects per treatment which are then normalize by n_samples."""
        self.avg_effects = {treatment:None for treatment in self.treatments}

        for treatment in self.treatments:
            #Grabbing all identified effects in a treatment
            df = pd.concat(self.snps_per_treatment[treatment])
            effects = set([effect.split(' ')[0] for effect in df['EFFECT'] \
                if not pd.isna(effect)])

            #Counts observed effect
            effects = dict(sorted({effect:0 for effect in effects}.items()))
            for effect in df['EFFECT']:
                if not pd.isna(effect):
                    effect = effect.split(' ')[0]
                    effects[effect] += 1
            
            #Normalizes counts by n_samples
            for effect,count in effects.items():
                effects[effect] = count/len(self.snps_per_treatment[treatment])

            self.avg_effects[treatment] = effects
    
    def get_gene_counts(self,dfs):
        """Returns the gene counts over subsetted samples. dfs are
        subsetter according to treatment, timepoint or microcosm.
        Results are normalized by n samples."""
        df = pd.concat(dfs).dropna()
        gene_counts = {gene:0 for gene in set(df['GENE'])}
        for gene in df['GENE']:
            gene_counts[gene] += 1
        gene_counts_per_sample = gene_counts
        for gene,count in gene_counts.items():
            gene_counts_per_sample[gene] = count/len(dfs)
        return gene_counts

    def get_gene_counts_per_treatment(self):
        """Gets the normalized gene count subsetted by treatment."""
        self.gene_counts_per_treatment = {treatment:None for treatment in self.treatments}
        for treatment in self.treatments:
            self.gene_counts_per_treatment[treatment] = \
                self.get_gene_counts(self.snps_per_treatment[treatment])

m = Mutations()
m.get_snps()
