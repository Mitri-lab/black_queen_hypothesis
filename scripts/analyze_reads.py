from Bio.SeqRecord import SeqRecord
from samples import Samples
from os.path import join,exists
from hgt import Hgt
import json
import pandas as pd
import pysam

s = Samples()


def get_origins():
    for strain, samples in s.strains.items():
        for sample in samples:
            if sample['platform'] == 'pacbio':
                print(sample['name'])
                hgt = Hgt(sample)
                hgt.chunk_assembly()
                hgt.mapper()
                hgt.get_mapping_stats()
                hgt.dump_origins()
                df = concat_origins(hgt.filtered)
                print(df)
def get_del(sample, chrom, pos):
    reads = []
    a = pysam.AlignmentFile(
        join(sample['dir_name'], 'aligned.sorted.bam'), 'rb')
    for read in a.fetch(chrom, pos, pos+1):
        if (not read.is_secondary) | (not read.is_unmapped) | (read.mapq == 60):
            reads.append(read)
    return reads

def concat_origins(origins):
    df = pd.DataFrame(columns=['chromosome','position','length','origins'])
    i = -1
    for name,contigs in origins.items():
        prev = 0
        prev_strains = None
        for pos, strains in contigs.items():
            pos = int(pos)
            if (pos -1 == prev) & (prev_strains == strains):
                df.at[i,'length'] += 1
            else:
                i += 1
                df.loc[i] = [name,pos,1,strains]
            prev = pos
            prev_strains = strains
    return df

def parse_origins():
    for strain, samples in s.strains.items():
        for sample in samples:
            if sample['platform'] == 'pacbio':
                f = join(sample['dir_name'],'origins.json')
                with open(f,'r') as handle:
                    data = handle.read()
                origins = json.loads(data)
                df = concat_origins(origins)
                for chrom,pos in zip(df['chromosome'],df['position']):
                    reads = get_del(sample,chrom,pos)
                    if len(reads) == 0:
                        print(sample['name'])
                        print(chrom,pos)
                        print(df)

def plot_insertions():
    samples = []
    for starin,samples in s.strains.items():
        for sample in samples:
            f = join(sample['dir_name'],'insertions.hgt.gbk')
            if exists(f):
                samples.append(sample)

