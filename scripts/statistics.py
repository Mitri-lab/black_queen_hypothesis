import pandas as pd
from scipy.stats import ttest_ind, kruskal
from diversity import get_mutations, get_variants
from snippy_core import distance_tree
from samples import Samples
import plotly.express as px
from os.path import join
s = Samples()


strains = {s.abbreviations['at']: 'At',
           s.abbreviations['ct']: 'Ct',
           s.abbreviations['oa']: 'Oa',
           s.abbreviations['ms']: 'Ms'}


def kruskal_treatment(df,column):
    out = pd.DataFrame(columns=['treatment1','treatment2','timepoint','t','p'])
    for j in ['T11', 'T22', 'T33', 'T44']:
        t2 = df[(df['timepoint'] == j) & (
            df['treatment'] == 2)][column].to_numpy()
        t3 = df[(df['timepoint'] == j) & (
            df['treatment'] == 3)][column].to_numpy()
        t4 = df[(df['timepoint'] == j) & (
            df['treatment'] == 4)][column].to_numpy()
        t, p = kruskal(t2, t3)
        out.loc[len(out)] = ['2','3',j,t,p]
        t, p = kruskal(t2, t4)
        out.loc[len(out)] = ['2','4',j,t,p]
    return out

def mutation_rate():
    df = get_mutations(add_T0=False)
    df = df[df['strain'] == strains[s.abbreviations['ct']]]
    fixation_rate = kruskal_treatment(df,'fixation_rate')
    fixation_rate = fixation_rate[fixation_rate['p'] <= 0.05]
    acc_rate = kruskal_treatment(df,'acc_rate')
    acc_rate = acc_rate[acc_rate['p'] <= 0.05]
    return fixation_rate,acc_rate

def mutations():
    df = get_mutations(add_T0=False)
    df = df[df['strain'] == strains[s.abbreviations['ct']]]
    fixed = kruskal_treatment(df,'fixed')
    fixed = fixed[fixed['p'] <= 0.05]
    total = kruskal_treatment(df,'mutations')
    total = total[total['p'] <= 0.05]
    ratio = kruskal_treatment(df,'fixed_total_ratio_number')
    ratio = ratio[ratio['p'] <= 0.05]
    f = join('..', 'variants', 'variants_comp_mapping.csv')
    df = get_variants(f,'illumina')
    df = df[df['strain'] == 'Ct']
    variants = kruskal_treatment(df,'hill')

    #variants = variants[variants['p'] <= 0.05]
    
    return total,fixed, variants,ratio

def genome_legnth():
    df = pd.read_csv(join('..','variants','assembly_length.csv'))
    df.insert(len(df.columns),'timepoint','T44')
    df = df[df['strain'] == 'Ct']
    out = kruskal_treatment(df,'deletions')

df = pd.read_csv(join('..','variants','genes.csv'))
df.insert(len(df.columns),'timepoint','T44')
df = df[df['strain'] == s.abbreviations['ct']]
out = kruskal_treatment(df,'genes')

def distance():
    dist = distance_tree('ct')
    df = pd.DataFrame(columns=['timepoint','treatment','distance'])
    for i,d in zip(dist.index[1:],dist['Ancestor'][1:]):
        tp = i[:3]
        t = int(i[5])
        df.loc[len(df)] = [tp,t,d]
    return kruskal_treatment(df,'distance')
#t,f,v,r = mutations()
#d = distance()
#l = genome_legnth()