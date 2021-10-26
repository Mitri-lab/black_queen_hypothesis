from samples import Samples
from samples import Experiment
import pandas as pd
from os.path import join
from os.path import exists
import glob

s = Samples()
e = Experiment()

class Stats():
    def get_sample_counts(self):
        stats = pd.DataFrame(0,columns=['treatment','T11','T22','T33','T44'],index=[1,2,3,4])
        for treatment in e.treatments:
            stats.at[treatment,'treatment'] = treatment
        names = []
        for strain,samples in s.strains.items():
            for sample in samples:
                names.append(sample['name'])
        counted = {name:False for name in set(names)}
        for strain,samples in s.strains.items():
            for sample in samples:
                if (not counted[sample['name']]) & (sample['platform'] == 'illumina'):
                    stats.at[sample['treatment'],sample['timepoint']] += 1
                    counted[sample['name']] = True
        stats.to_csv('sample_counts.csv',index=False)

    def get_read_mappings(self):
        for timepoint in e.timepoints:
            for treatment in [3,4]:
                names = []
                for strain,samples in s.strains.items():
                    for sample in samples:
                        if (sample['platform'] == 'illumina') & (sample['treatment'] == treatment) \
                            & (sample['timepoint'] == timepoint):
                            names.append(sample['name'])
                names = sorted(set(names))
                if treatment == 3:
                    out = pd.DataFrame(columns=names,index=['ct','at','ms'])
                if treatment == 4:
                    out = pd.DataFrame(columns=names,index=['ct','at','ms','oa'])
                for strain,samples in s.strains.items():
                    for sample in samples:
                        if (sample['name'] in names) & (sample['platform'] == 'illumina'):
                            stats = join(sample['dir_name'],'flagstat.tsv')
                            df = pd.read_csv(stats,sep='\t',names=['stats','unknown','name'])
                            df.index = df['name'].to_list()
                            out.at[s.abbreviations[sample['strain']],sample['name']] = df.loc['mapped %']['stats']
                df.index.name = 'strain'
                out.to_csv(timepoint+'_'+str(treatment)+'.csv')

    def get_ct_snps(self):
        out = pd.DataFrame(0,columns=['treatment 2','treatment 3','treatment 4'],index=['ct'])
        for sample in s.strains[s.abbreviations['ct']]:
            if sample['platform'] == 'illumina':
                f = join(sample['dir_name'],'snippy','snps.tab')
                df = pd.read_csv(f,sep='\t')
                out.at['ct','treatment '+str(sample['treatment'])] += len(df)
        out.to_csv('at_snps.csv',index=False)

    def get_ct_effects(self):
        work = '/work/FAC/FBM/DMF/smitri/evomicrocomm/genome_size/data/'
        effects = pd.concat([pd.read_csv(csv,sep='\t') for csv \
            in glob.glob(join(work,'T*','ct','snippy','snps.tab'))]).dropna()['EFFECT']
        effects = set([effect.split(' ')[0] for effect in effects])
        out = pd.DataFrame(0,columns=['treatment 2','treatment 3','treatment 4'],index=effects)
        for sample in s.strains[s.abbreviations['ct']]:
            if sample['platform'] == 'illumina':
                f = join(sample['dir_name'],'snippy','snps.tab')
                if exists(f):
                    effects = pd.read_csv(f,sep='\t').dropna()['EFFECT']
                    for effect in effects:
                        out.at[effect.split(' ')[0],'treatment '+str(sample['treatment'])] += 1
        out.index.name = 'effect'
        out.to_csv('ct_effects.csv')
    
    def get_crp(self):
        work = '/work/FAC/FBM/DMF/smitri/evomicrocomm/genome_size/data/'
        dfs = []
        for sample in s.strains[s.abbreviations['at']]:
            if sample['platform'] == 'illumina':
                df = pd.read_csv(join(sample['dir_name'],'snippy','snps.tab'),sep='\t')
                if len(df) > 0:
                    df.index = [sample['name']]*len(df)
                    dfs.append(df)
        out = pd.concat(dfs)
        out = out[out['GENE'] == 'crp']
        out = out[['CHROM','POS','TYPE','EFFECT','GENE']]
        out['EFFECT'] = [effect.split(' ')[0] for effect in out['EFFECT']]
        out.index.name = 'sample'
        print(out)
        out.to_csv('crp.csv')

    def get_no_alignments(self):
        dfs = []
        for sample in s.strains[s.abbreviations['at']]:
            if sample['platform'] == 'illumina':
                df = pd.read_csv(join(sample['dir_name'],'no_alignment_regions.tsv'),sep='\t').dropna()
                df = df[df['product'] != 'hypothetical protein']
                if len(df) > 0:
                    df.index = [sample['name']] * len(df)
                    dfs.append(df)
        return pd.concat(dfs)



stats = Stats()
o = stats.get_no_alignments()
#m = Mutations()
#m.get_snps()
