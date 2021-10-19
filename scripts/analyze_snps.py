from samples import Samples
import pandas as pd
from os.path import join
from os.path import exists

s = Samples()

class Snps():
    def avg_snps_per_treatment(self,strain):
        self.avg_snps = dict()
        snps = dict()
        for sample in s.strains[strain]:
            if (sample['platform'] == 'illumina') & exists(join(sample['dir'],'snippy','snps.tab')):
                snps[sample['treatment']] = []
                self.avg_snps[sample['treatment']] = None
        
        for sample in s.strains[strain]:
            if (sample['platform'] == 'illumina') & exists(join(sample['dir'],'snippy','snps.tab')):
                df = pd.read_csv(join(sample['dir'],'snippy','snps.tab'),sep='\t')
                snps[sample['treatment']].append(len(df))

        for treatment in snps.keys():
            self.avg_snps[treatment] = sum(snps[treatment])/len(snps[treatment])

for strain in s.strains:
    p = Snps()
    p.get_snps_per_treatment(strain)
    s.strains[strain] = p