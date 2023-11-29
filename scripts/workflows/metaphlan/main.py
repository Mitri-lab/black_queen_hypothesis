from os import mkdir
from os.path import join
from os.path import exists
from subprocess import call
from samples import Samples

s = Samples()
work = '/work/FAC/FBM/DMF/smitri/evomicrocomm/genome_size/data/'

def call_metaphlan():
    names = []
    for strain,samples in s.strains.items():
        for sample in samples:
            if sample['platform'] == 'illumina':
                names.append(sample['name'])

    for name in set(names):
        f = join(work,name,name+'.pkl')
        if not exists(f):
            call(['sbatch','metaphlan.sh',name])

def create_db_markers():
    for strain in s.strains:
        fname = 's__' + strain.replace(' ','_') + '.fna'
        if not exists(fname):
            call(['bash','db_markers.sh',strain.replace(' ','_'),s.abbreviations[strain]])

def strainphlan_timepoint(timepoint):
    for strain,samples in s.strains.items():
        d = join(work,s.abbreviations[strain],'strainphlan')
        if not exists(d):
            mkdir(d)
        regex = timepoint + '*'
        call(['sbatch','strainphlan.sh',regex,strain,s.abbreviations[strain]])