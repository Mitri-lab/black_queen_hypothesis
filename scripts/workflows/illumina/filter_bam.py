import pysam
from os.path import join, split
import sys
from samples import Samples
from os import symlink

s = Samples()


def get_strains(d):
    strains = []
    for strain, samples in s.strains.items():
        for sample in samples:
            if sample['name'] == d:
                strains.append(s.abbreviations[strain])
    return strains


def symlink_bam(b):
    t = join(split(b)[0], 'mapped_reads.filtered.sorted.bam')
    symlink(b, t)


def filter_bam(d, b):
    q_strain = split(split(b)[0])[1]
    ref_strains = [i for i in get_strains(d) if i != q_strain]
    print(ref_strains)
    if len(ref_strains) == 0:
        symlink_bam(b)
        return
    o = join(split(b)[0], 'mapped_reads.filtered.sorted.bam')

    ids = []
    for j in ref_strains:
        print(j)
        f = join(s.work, d, j, 'read_ids.txt')
        with open(f, 'r') as handle:
            data = list(set(handle.read().split(',')))
            ids += data
    id_dict = {key: True for key in set(ids)}
    a = pysam.AlignmentFile(b)
    filtered = pysam.AlignmentFile(o, 'wb', template=a)
    for read in a:
        try:
            if id_dict[read.qname]:
                pass
        except KeyError:
            filtered.write(read)


if __name__ == "__main__":
    d = sys.argv[1]
    b = sys.argv[2]
    print(d,b)
    filter_bam(d, b)
