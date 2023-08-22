from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
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
s = Samples()
colors = {'ct': ['#c0369d', '#fa7876', '#6a3f99'],
          'at': ['#1f9e75', '#fa7876', '#6a3f99'],
          'all': ['#fa7876', '#c0369d', '#6a3f99']}


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
        if node.name[5] == '1':
            name_face.background.color = '#a493c7'
            name_face.border.color = '#4b2991'
        if node.name[5] == '2':
            name_face.background.color = '#df9acd'
            name_face.border.color = '#c0369d'
        if node.name[5] == '3':
            name_face.background.color = '#fcbbba'
            name_face.border.color = '#c0369d'
        if node.name[5] == '4':
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
    cmd = ['run_gubbins.py', '-p', 'gubbins', '--outgroup',
           'Reference', prefix+'.clean.full.aln']
    call(' '.join(cmd), shell=True)
    if exists('gubbins.final_tree.tre'):
        style_tree(specie, prefix)
    else:
        print('Tree build failed.')
    os.chdir(join('/', 'users', 'eulrich', 'black_queen_hypothesis', 'scripts'))


#caller('all', ['T44'], [2, 3, 4], 'ct', 'illumina', 'ct')
def font_size(fig):
    """Style function for figures setting fot size and true black color."""
    j = 7
    fig.update_layout(font={'size': j, 'color': 'black'})
    for a in fig['layout']['annotations']:
        a['font']['size'] = j
        a['font']['color'] = 'black'
    fig['layout']['title']['font']['size'] = j
    fig['layout']['title']['font']['color'] = 'black'
    fig['layout']['legend']['title']['font']['size'] = j
    fig['layout']['legend']['title']['font']['color'] = 'black'
    fig.for_each_xaxis(lambda axis: axis.title.update(
        font=dict(size=j, color='black')))
    fig.for_each_yaxis(lambda axis: axis.title.update(
        font=dict(size=j, color='black')))
    fig.update_layout(
        margin=dict(l=55, r=60, t=0, b=20,autoexpand=False),coloraxis_colorbar=dict(len=0.6,thickness=10))

    return fig

def distance_tree(abb):
    w = 300
    h = 300
    f = join(s.work, abb, 'trees', 'gubbins.filtered_polymorphic_sites.fasta')
    #f = join(s.work, abb, 'trees', abb+'.aln')
    aln = AlignIO.read(open(f), 'fasta')
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(aln)
    s.labels['Reference'] = 'Ancestor'
    labels = [s.labels[i] for i in dm.names]
    dm.names = labels
    df = pd.DataFrame(columns=dm.names)
    for j, i in enumerate(dm):
        df.loc[j] = i
    df.index = dm.names
    fig = px.imshow(df,width=w,height=h)
    fig = font_size(fig)
    fig.write_image(join('..', 'plots', 'plots', 'distance_matrix.svg'))
    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(dm)
    tree.root_with_outgroup('Ancestor')
    f_out = join(s.work, abb, 'trees', 'identity.tre')
    with open(f_out, 'w') as handle:
        handle.write(tree.format('newick'))
    style_tree(abb, abb, f_out)



# caller('all', ['T11','T22', 'T33', 'T44'], [1, 3, 4], 'ct', 'illumina', 'at')
