from Bio.SeqRecord import SeqRecord
from samples import Samples
from os.path import join, exists
from Bio import SeqIO
from hgt import Hgt
import json
import pandas as pd
import pysam
from dna_features_viewer import GraphicFeature, GraphicRecord

s = Samples()


def get_origins():
    """Calls hgt class on every sample"""
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

def get_origins_anceteral():
    strains = ['at','ct','oa','ms']
    for strain in strains:
        sample = dict()
        sample['dir_name'] = join(s.work,strain)
        sample['name'] = strain
        sample['strain'] = s.abbreviations[strain]
        hgt = Hgt(sample)
        hgt.chunk_assembly()
        hgt.mapper()
        hgt.get_mapping_stats()
        hgt.dump_origins()
        df = concat_origins(hgt.filtered)
        print(df)

def get_del(sample, chrom, pos):
    """For insertion detection ancesteral strain is already
    mapped to the assembly of the evolved sample.
    We can chech the foreing origins and see if anceteral sequences
    are present at evolved assembly position. If yes, foreign sequence
    was already present in anceteral strain. If not, it's likeyly
    that a HGT event was detected taking place during evolution
    experiment."""
    reads = []
    a = pysam.AlignmentFile(
        join(sample['dir_name'], 'aligned.sorted.bam'), 'rb')
    for read in a.fetch(chrom, pos, pos+1):
        if (not read.is_secondary) | (not read.is_unmapped) | (read.mapq == 60):
            reads.append(read)
    return reads


def concat_origins(origins):
    """Concats origins. In input every position has strain info.
    This function concatenates positions following eachother and
    outputs the start position and the length of the foreing origin."""
    df = pd.DataFrame(columns=['chromosome', 'position', 'length', 'origins'])
    i = -1
    for name, contigs in origins.items():
        # Previous pos
        prev = 0
        # Previous strain
        prev_strains = None
        for pos, strains in contigs.items():
            pos = int(pos)
            if (pos - 1 == prev) & (prev_strains == strains):
                df.at[i, 'length'] += 1
            else:
                i += 1
                df.loc[i] = [name, pos, 1, strains]
            prev = pos
            prev_strains = strains
    return df


def parse_origins():
    """Reads dumped jsons"""
    for strain, samples in s.strains.items():
        for sample in samples:
            if sample['platform'] == 'pacbio':
                f = join(sample['dir_name'], 'origins.json')
                with open(f, 'r') as handle:
                    data = handle.read()
                origins = json.loads(data)
                df = concat_origins(origins)
                for chrom, pos in zip(df['chromosome'], df['position']):
                    reads = get_del(sample, chrom, pos)
                    if len(reads) == 0:
                        print(sample['name'])
                        print(chrom, pos)
                        print(df)


def plot_insertions():
    """Plots foreing origins"""
    for strain, samples in s.strains.items():
        for sample in samples:
            if sample['platform'] == 'pacbio':
                f = join(sample['dir_name'], 'insertions.hgt.gbk')
                if exists(f):
                    features = []
                    contig = [contig for contig in SeqIO.parse(f, 'genbank')][0]
                    for feature in contig.features:
                        if feature.type == 'CDS':
                            gf = GraphicFeature(start=feature.location.start, end=feature.location.end, strand=feature.location.strand,\
                                    color="#ffd700",label=feature.qualifiers['product'][0])
                            features.append(gf)
                    record = GraphicRecord(sequence_length=len(contig),features=features)
                    record.plot_on_multiple_pages(join(sample['dir_name'],
                    "multipage_plot.pdf"),
                    nucl_per_line=5000,
                    lines_per_page=10,
                    plot_sequence=False
                )
                    
plot_insertions()