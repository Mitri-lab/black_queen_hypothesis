from Bio.Phylo.TreeConstruction import DistanceCalculator,DistanceTreeConstructor
from Bio import AlignIO
import plotly.express as px
import dendropy
import numpy as np
from Bio import Phylo
import itertools
from samples import Samples
import pandas as pd
from os.path import join, exists, split
from os import symlink, remove, unlink, mkdir, chdir
from shutil import move, copyfile
from subprocess import call, DEVNULL, STDOUT
from glob import glob
from ete4 import Tree, TreeStyle, AttrFace, faces
import os
from Bio import SeqIO
os.environ['QT_QPA_PLATFORM'] = 'offscreen'
import itertools
s = Samples()


def leave_names():
    for specie, samples in s.strains.items():
        for sample in samples:
            if (sample['platform'] == 'illumina') | (sample['platform'] == 'pacbio'):
                src = join(sample['dir_name'], 'snippy')
                trgt = join(sample['dir_name'], sample['name'])
                if exists(trgt):
                    unlink(trgt)
                    symlink(src, trgt, target_is_directory=True)
                else:
                    unlink(trgt)
                    symlink(src, trgt, target_is_directory=True)


def backup_aligned():
    for specie, samples in s.strains.items():
        for sample in samples:
            if (sample['platform'] == 'illumina') | (sample['platform'] == 'pacbio'):
                f = join(sample['dir_name'], 'snippy', 'snps.aligned.fa')
                d = join(s.work, 'bulk', sample['name'])
                if not exists(d):
                    mkdir(d)
                d = join(d, s.abbreviations[sample['strain']])
                if not exists(d):
                    mkdir(d)
                t = join(d, sample['name']+'.fa')
                if not exists(t):
                    copyfile(f, t)


def restore_aligned():
    for specie, samples in s.strains.items():
        for sample in samples:
            if (sample['platform'] == 'illumina'):
                src = join(s.work, 'bulk',
                           sample['name'], sample['name'] + '.fa')
                t = join(sample['dir_name'], 'snippy', 'snps.aligned.fa')
                if exists(t):
                    remove(t)
                    print(src, t)
                    copyfile(src, t)


def filter_aligned():
    for specie, samples in s.strains.items():
        for sample in samples:
            if (sample['platform'] == 'illumina'):
                out = []
                f = join(sample['dir_name'], 'snippy', 'snps.aligned.fa')
                contigs = [contig for contig in SeqIO.parse(f, 'fasta')]
                for c in contigs:
                    if s.abbreviations[sample['strain']] in c.name:
                        c.name, c.description = '', ''
                        out.append(c)
                with open(f, 'w') as handle:
                    SeqIO.write(out, handle, 'fasta')


def my_layout(node):
    if node.is_leaf():
        name_face = AttrFace(
            "name", fsize=10, ftype='Open Sans', fgcolor='black')
        faces.add_face_to_node(name_face, node, column=0,
                               position="branch-right")
        if node.name[4] == '1':
            name_face.background.color = '#a493c7'
            name_face.border.color = '#4b2991'
        if node.name[4] == '2':
            name_face.background.color = '#df9acd'
            name_face.border.color = '#c0369d'
        if node.name[4] == '3':
            name_face.background.color = '#fcbbba'
            name_face.border.color = '#c0369d'
        if node.name[4] == '4':
            name_face.background.color = '#f5ebd0'
            name_face.border.color = '#edd9a3'
        name_face.border.width = 1


def style_tree(specie, prefix,f=None):
    if not f:
        f = join(s.work, specie, 'trees', 'gubbins.node_labelled.final_tree.tre')
    tree = Tree(f, format=1)
    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.layout_fn = my_layout
    #ts.branch_vertical_margin = -1
    #ts.scale = 20
    out = join(s.work, specie, 'trees', prefix+'.svg')
    tree.render(out, tree_style=ts)


def caller(cosm, timepoint, treatment, specie, platform, prefix):
    print(prefix)
    os.chdir(join(s.work, specie, 'trees'))
    ref = join(s.work, 'ancestors', specie, 'bakta', 'snippy.gbff')
    cmd = ['snippy-core', '--ref', ref, '--prefix', prefix]
    for sample in s.strains[s.abbreviations[specie]]:
        if sample['platform'] == platform:
            if sample['treatment'] in treatment:
                if cosm == 'all':
                    if timepoint == 'all':
                        cmd.append(join(sample['dir_name'], sample['name']))
                    else:
                        if sample['timepoint'] in timepoint:
                            cmd.append(
                                join(sample['dir_name'], sample['name']))
                else:
                    if sample['cosm'] == cosm:
                        if timepoint == 'all':
                            cmd.append(
                                join(sample['dir_name'], sample['name']))
                        else:
                            if sample['timepoint'] == timepoint:
                                cmd.append(
                                    join(sample['dir_name'], sample['name']))
    call(' '.join(cmd), shell=True)
    cmd = ['snippy-clean_full_aln', prefix +
           '.full.aln', '>', prefix+'.clean.full.aln']
    call(' '.join(cmd), shell=True)
    cmd = ['run_gubbins.py', '-p', 'gubbins','--outgroup', 'Reference', prefix+'.aln']
    call(' '.join(cmd), shell=True)
    if exists('gubbins.final_tree.tre'):
        style_tree(specie, prefix)
    else:
        print('Tree build failed.')
    os.chdir(join('/', 'users', 'eulrich', 'black_queen_hypothesis', 'scripts'))


def newick_to_matrix(file):
    t = Phylo.read(file, 'newick')

    d = {}
    for x, y in itertools.combinations(t.get_terminals(), 2):
        v = t.distance(x, y)
        d[x.name] = d.get(x.name, {})
        d[x.name][y.name] = v
        d[y.name] = d.get(y.name, {})
        d[y.name][x.name] = v
    for x in t.get_terminals():
        d[x.name][x.name] = 0

    m = pd.DataFrame(d)
    return m


#caller('all', ['T44'], [2, 3, 4], 'ct', 'illumina', 'ct')

def distance_tree(f):
    aln = AlignIO.read(open(f), 'fasta')
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(aln)
    df = pd.DataFrame(columns=dm.names)
    for j, i in enumerate(dm):
        df.loc[j] = i
    df.index = dm.names
    fig = px.imshow(df)
    fig.write_image(join('..','plots','plots','distance_matrix.svg'))
    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(dm)
    tree.root_with_outgroup('Reference')
    f_out = join(s.work, 'ct', 'trees', 'identity.tre')
    with open(f_out,'w') as handle:
        handle.write(tree.format('newick'))
    style_tree('ct','ct',f_out)
    return df

df = distance_tree(f = join(s.work, 'ct', 'trees', 'ct.aln'))

"""combs = itertools.combinations(df.columns.to_list(),2)
out  = pd.DataFrame(columns=['dist','treatment'])
samples = list(df.columns)
samples.remove('Reference')
combs = itertools.combinations(samples,2)
dist = {'2':[],'3':[],'4':[]}
for (p1,p2) in combs:
    if p1[4] != '2':
        break
    else:
        d = tree.distance(p1,p2)
        dist[p2[4]].append(d)
        out.loc[len(out)] = [d,p2[4]]

g1 = ['T44.2.1','T44.2.2.re','T44.2.3.re','T44.2.5.re']
g2 = ['T44.3.1','T44.2.3','T44.3.3.re','T44.3.4']"""

# caller('all', ['T11','T22', 'T33', 'T44'], [1, 3, 4], 'ct', 'illumina', 'at')

