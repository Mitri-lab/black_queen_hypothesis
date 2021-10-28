from samples import Samples
from samples import Experiment
from plotting import Plotting
import pandas as pd
from os.path import join
from os.path import exists

s = Samples()
e = Experiment()
p = Plotting()

def plot_genes():
    """This plots the mutated genes in n microocms."""
    for strain,samples in s.strains.items():
        treatments = s.treatments[strain]
        samples = [sample for sample in samples if sample['platform'] == 'illumina']
        for treatment in treatments:
            genes = []
            for sample in samples:
                if sample['treatment'] == treatment:
                    snps = join(sample['dir_name'],'snippy','snps.tab')
                    if exists(snps):
                        genes += pd.read_csv(snps,sep='\t').dropna()['GENE'].to_list()
            out = pd.DataFrame(0,columns=e.timepoints,index=set(genes))
            for sample in samples:
                if sample['treatment'] == treatment:
                    snps = join(sample['dir_name'],'snippy','snps.tab')
                    if exists(snps):
                        for gene in set(pd.read_csv(snps,sep='\t').dropna()['GENE']):
                            out.at[gene,sample['timepoint']] += 1 
            fig = p.trajectories(out)
            title = [strain,'in','treatment',str(treatment)]
            fig.update_layout(
                    title = ' '.join(title),
                    xaxis_title = 'timepoint',
                    yaxis_title = 'observed in n microcosms'
                )
            fig.write_image(join('..','plots','genes',' '.join(title).replace(' ','_')+'.png'))

plot_genes()