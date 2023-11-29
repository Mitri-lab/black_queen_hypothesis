from samples import Samples
from subprocess import call
from os.path import join as d_join, exists
import json
s = Samples()


def caller():
    samples = []
    timepoint = 'T44'
    species = [s.abbreviations[specie] for specie in ['at']]
    for specie in species:
        specie = s.abbreviations[specie]
        samples.append(d_join(s.work, 'ancestrals', 'ct', 'reads.sig'))

    for specie in species:
        for sample in s.strains[specie]:
            # & (sample['timepoint'] == timepoint):
            if (sample['platform'] == 'illumina'):
                f = d_join(sample['dir_name'], 'reads.sig')
                if exists(f):
                    samples.append(f)
    out = d_join(s.work, 'bulk', 'sourmash')
    samples = sorted(samples, reverse=True)
    print(samples)
    call(' '.join(['sourmash', 'compare', '-p', '31',
         ' '.join(samples), '-o', out]), shell=True)


def cosm_caller():
    samples = []
    # samples.append(d_join(s.work, 'ancestors', specie, 'reads.sig'))
    for specie in ['at', 'ct']:
        for sample in s.strains[s.abbreviations[specie]]:
            # & (sample['timepoint'] == 'T44'):
            if (sample['platform'] == 'illumina') & (sample['timepoint'] == 'T44'):
                f = d_join(sample['dir_name'], 'reads.sig')
                samples.append(f)
    out = d_join(s.work, specie, 'sourmash')
    samples = sorted(samples, reverse=True)
    print(samples)
    call(' '.join(['sourmash', 'compare', '-p', '16', '-k', '31',
         ' '.join(samples), '-o', out, '--csv', out+'.csv']), shell=True)


def genome_content():
    pass


def plot(species):
    out = d_join(s.work, species)
    f = d_join(s.work, species, 'sourmash')
    call(' '.join(['sourmash', 'plot','--pdf',
         '--labels', '--output-dir', out, f]), shell=True)


def fix_labels(species):
    timepoint = 'T44'
    for sample in s.strains[s.abbreviations[species]]:
        if (sample['platform'] == 'illumina'):
            f = d_join(sample['dir_name'], 'reads.sig')
            if exists(f):
                print(f)
                f = open(f)
                data = json.load(f)
                data[0]['name'] = '_'.join(
                    [species[0].upper() + species[1], str(sample['treatment']), str(sample['cosm'])])
                with open(d_join(sample['dir_name'], 'reads.sig'), 'w') as out:
                    json.dump(data, out)
            else:
                pass


def pacbio(specie):
    samples = []
    samples.append(d_join(s.work, 'ancestors', specie,
                          'assembly.contigs.polypolish.sig'))

    for sample in s.strains[s.abbreviations[specie]]:
        if (sample['platform'] == 'pacbio'):
            f = d_join(sample['dir_name'], 'assembly.sig')
            samples.append(f)
    out = d_join(s.work, specie, 'sourmash_pacbio')
    samples = sorted(samples, reverse=True)
    print(samples)
    call(' '.join(['sourmash', 'compare', '-p', '16',
         ' '.join(samples), '-o', out]), shell=True)


def fix_labels_pacbio(species):
    for sample in s.strains[s.abbreviations[species]]:
        if (sample['platform'] == 'pacbio'):
            f = open(d_join(sample['dir_name'], 'assembly.sig'))
            data = json.load(f)
            data[0]['name'] = '_'.join(
                [species[0].upper() + species[1], sample['treatment'], sample['cosm']])
            '_' + species[0].upper() + species[1]
            with open(d_join(sample['dir_name'], 'assembly.sig'), 'w') as out:
                json.dump(data, out)


def plot_pacbio(specie):
    out = d_join(s.work, specie)
    f = d_join(s.work, specie, 'sourmash_pacbio')
    call(' '.join(['sourmash', 'plot', '--pdf',
         '--labels', '--output-dir', out, f]), shell=True)


def pacbio_caller():
    for sp in ['at', 'ct', 'ms', 'oa']:
        fix_labels_pacbio(sp)
        pacbio(sp)
        plot_pacbio(sp)


"""q2 = d_join(s.work, 'T44.4.4', 'ct', 'reads.sig')
r = d_join(s.work, 'ancestors', 'ct', 'assembly.contigs.polypolish.sig')
q = d_join(s.work, 'T44.2.1', 'ct', 'reads.sig')
cmd = ['sourmash', 'search', '--ksize', '31', q, q2, '--containment']
call(' '.join(cmd), shell=True)
"""
