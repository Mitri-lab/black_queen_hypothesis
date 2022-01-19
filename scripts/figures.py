from dna_features_viewer import GraphicFeature, GraphicRecord
from samples import Samples
import pysam
from os.path import join

s = Samples()
sample = s.strains['Agrobacterium tumefaciens'][9]


def plot_chunks():
    reads = []
    c = 'ctg.s2.34.arrow'
    step = 100
    window = 500
    step_ids = (586, 636)
    n_steps = step_ids[1] - step_ids[0]
    steps = range(step_ids[0], step_ids[1])
    length = step_ids[1] * step + window
    features = []
    starts = []
    ends = []
    for i in steps:
        
        label = '.'.join([c, str(i)])
        reads.append(label)
        if i*step+window > length:
            break
        else:
            gf = GraphicFeature(start=i*step, end=i*step+window, strand=1,
                                color="#ffd700", label=label)
            starts.append(i*step)
            ends.append(i*step+window)
            features.append(gf)
    print(sorted(ends)[-1]-sorted(starts)[0])
    record = GraphicRecord(sequence_length=sorted(ends)[-1], features=features)
    record = record.crop((sorted(starts)[0],sorted(ends)[-1]))
    record.plot_on_multiple_pages("chunks.pdf",
                                  nucl_per_line=5400,
                                  lines_per_page=10,
                                  plot_sequence=False
                                  )
    return reads

def plot_alignments():
    reads = plot_chunks()
    strains = ['at','oa']
    for strain in strains:
        a_reads = []
        features = []

        a = pysam.AlignmentFile(join(sample['dir_name'], 'hgt', strain+'.sam'), 'rb')
        for read in a:
            if read.qname in reads:
                if (not read.is_unmapped) & (not read.is_secondary) & (read.mapq == 60):
                    a_reads.append(read)
        starts = []
        ends = []
        for read in a_reads:
            if read.is_reverse:
                strand = -1
            else:
                strand = 1
            gf = GraphicFeature(start=read.reference_start, end=read.reference_end, strand=strand,
                                color="#ffd700", label=read.qname)
            starts.append(read.reference_start)
            ends.append(read.reference_end)
            features.append(gf)
        record = GraphicRecord(sequence_length=sorted(ends)[-1]+5000, features=features)
        record = record.crop((sorted(starts)[0]-2000,sorted(ends)[-1]+2500))
        record.plot_on_multiple_pages("oa.pdf",
                                    nucl_per_line=5000,
                                    lines_per_page=10,
                                    plot_sequence=False
                                    )


def plot_oa():
    reads = plot_chunks()
    a_reads = []
    features = []

    a = pysam.AlignmentFile(join(sample['dir_name'], 'hgt', 'oa.sam'), 'rb')
    for read in a:
        if read.qname in reads:
            if (not read.is_unmapped) & (not read.is_secondary) & (read.mapq == 60):
                a_reads.append(read)
    starts = []
    ends = []
    for read in a_reads:
        if read.is_reverse:
            strand = -1
        else:
            strand = 1
        gf = GraphicFeature(start=read.reference_start, end=read.reference_end, strand=strand,
                            color="#ffd700", label=read.qname)
        starts.append(read.reference_start)
        ends.append(read.reference_end)
        features.append(gf)
    record = GraphicRecord(sequence_length=sorted(ends)[-1]+5000, features=features)
    record = record.crop((sorted(starts)[0]-2000,sorted(ends)[-1]+2500))
    record.plot_on_multiple_pages("oa.pdf",
                                nucl_per_line=5000,
                                lines_per_page=10,
                                plot_sequence=False
                                )

def plot_at():
    reads = plot_chunks()
    a_reads = []
    features = []

    a = pysam.AlignmentFile(join(sample['dir_name'], 'hgt', 'at.sam'), 'rb')
    for read in a:
        if read.qname in reads:
            if (not read.is_unmapped) & (not read.is_secondary) & (read.mapq == 60):
                a_reads.append(read)
    starts = []
    ends = []
    for read in a_reads:
        if read.is_reverse:
            strand = -1
        else:
            strand = 1
        gf = GraphicFeature(start=read.reference_start, end=read.reference_end, strand=strand,
                            color="#ffd700", label=read.qname)
        starts.append(read.reference_start)
        ends.append(read.reference_end)
        features.append(gf)
    record = GraphicRecord(sequence_length=sorted(ends)[-1], features=features)
    record = record.crop((sorted(starts)[0],sorted(ends)[-1]))
    record.plot_on_multiple_pages("at.pdf",
                                nucl_per_line=5400,
                                lines_per_page=10,
                                plot_sequence=False
                                )

