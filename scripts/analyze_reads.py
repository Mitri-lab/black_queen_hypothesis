from re import I
from Bio.SeqRecord import SeqRecord
from samples import Samples
from os.path import join
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

def get_del(sample, chrom, pos):
    reads = []
    a = pysam.AlignmentFile(
        join(sample['dir_name'], 'aligned.sorted.bam'), 'rb')
    for read in a.fetch(chrom, pos, pos+1):
        reads.append(read)
    return reads

for strain, samples in s.strains.items():
    for sample in samples:
        if sample['platform'] == 'pacbio':
            f = join(sample['dir_name'],'origins.json')
            with open(f,'r') as handle:
                data = handle.read()
            origins = json.loads(data)
            df = pd.DataFrame(columns=['chromosome','position','length','origins'])
            i = -1
            prev = 0
            for name,contigs in origins.items():
                for pos, strains in contigs.items():
                    pos = int(pos)
                    if pos -1 == prev:
                        df.at[i,'chromosome'] = name
                        df.at[i,'position'] = pos
                        df.at[i,'length'] += 1
                        df.at[i,'origins'] = strains
                        prev = pos
                    else:
                        i += 1
                        df.at[i,'chromosome'] = name
                        df.at[i,'position'] = pos
                        df.at[i,'length'] = 1
                        df.at[i,'origins'] = strains
                        prev = pos 
                    pass
            print(sample['name'])
            print(df)