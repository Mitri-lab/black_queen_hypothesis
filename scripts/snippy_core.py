from samples import Samples
import pandas as pd
from os.path import join, exists, split
from os import symlink, remove, unlink, mkdir
from shutil import move, copyfile
from subprocess import call, DEVNULL, STDOUT
from glob import glob
from ete4 import Tree, TreeStyle, AttrFace, faces
import os
from Bio import SeqIO
os.environ['QT_QPA_PLATFORM'] = 'offscreen'

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
        


def style_tree(specie,prefix):
    f = join(s.work,specie,'trees','gubbins.node_labelled.final_tree.tre')
    tree = Tree(f, format=1)
    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.layout_fn = my_layout
    ts.branch_vertical_margin = -1
    ts.scale = 20
    out = join(s.work,specie,'trees',prefix+'.svg')
    tree.render(out, tree_style=ts)

def caller(cosm, timepoint, treatment, specie, platform, prefix):
    print(prefix)
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
    print(cmd)
    call(' '.join(cmd), shell=True)
    cmd = ['snippy-clean_full_aln', prefix +
           '.full.aln', '>', prefix+'.clean.full.aln']
    call(' '.join(cmd), shell=True, stdout=DEVNULL, stderr=STDOUT)
    cmd = ['run_gubbins.py', '-p', 'gubbins',
           '--outgroup', 'Reference', prefix+'.clean.full.aln']
    call(' '.join(cmd), shell=True, stdout=DEVNULL, stderr=STDOUT)
    if exists('gubbins.final_tree.tre'):
        style_tree(specie,prefix)
    else:
        print('Tree build failed.')
    fs = glob(prefix+'*') + glob('gubbins*')
    t = join(s.work, specie, 'trees')
    for f in fs:
        move(f, join(t, split(f)[-1]))


#caller('all', ['T44'], [2,3,4], 'ct', 'illumina', 'ct')
#caller('all', ['T11','T22', 'T33', 'T44'], [2, 3, 4], 'ct', 'illumina', 'ct')
