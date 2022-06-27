from samples import Samples
from subprocess import call
from os.path import join as d_join,exists
import json
s = Samples()


def caller():
    samples = []
    timepoint = 'T44'
    species = [s.abbreviations[specie] for specie in ['at','ms']]
    for specie in species:
        specie = s.abbreviations[specie]
        samples.append(d_join(s.work,specie,specie+'.sig'))

    for specie in species:
        for sample in s.strains[specie]:
            if (sample['platform'] == 'illumina') & (sample['timepoint'] == timepoint):
                f = d_join(sample['dir_name'],'reads.sig')
                samples.append(f)
    out = d_join(s.work,'bulk','sourmash')
    samples = sorted(samples,reverse=True)
    print(samples)
    call(' '.join(['sourmash','compare','-p','16',' '.join(samples),'-o',out]),shell=True)

def cosm_caller():
    samples = []
    treatment = 4
    cosm = 3
    species = [s.abbreviations[specie] for specie in ['at']]
    for specie in species:
        specie = s.abbreviations[specie]
        samples.append(d_join(s.work,specie,specie+'.sig'))

    for specie in species:
        for sample in s.strains[specie]:
            if (sample['platform'] == 'illumina') & (sample['treatment'] == treatment) & (sample['cosm'] == cosm):
                f = d_join(sample['dir_name'],'reads.sig')
                samples.append(f)
    out = d_join(s.work,s.abbreviations[specie],'sourmash')
    samples = sorted(samples,reverse=True)
    print(samples)
    call(' '.join(['sourmash','compare','-p','16',' '.join(samples),'-o',out]),shell=True)

def plot(species):
    out = d_join(s.work,species)
    f = d_join(s.work,species,'sourmash')
    call(' '.join(['sourmash','plot','--pdf','--labels','--output-dir',out,f]),shell=True)

def fix_labels(species):
    timepoint = 'T44'
    for sample in s.strains[s.abbreviations[species]]:
        if (sample['platform'] == 'illumina'):
            f = open(d_join(sample['dir_name'],'reads.sig'))
            data = json.load(f)
            data[0]['filename'] = sample['name'] + '_' + species[0].upper() + species[1]
            with open(d_join(sample['dir_name'],'reads.sig'),'w') as out:
                json.dump(data,out)
