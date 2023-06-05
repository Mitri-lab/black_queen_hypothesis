import pysam
from os.path import join,split
import sys

def dump_read_ids(f):
    reads = []
    a = pysam.AlignmentFile(f)
    for read in a:
        reads.append(read.qname)
    o = join(split(f)[0],'read_ids.txt')
    with open(o, 'w') as handle:
        handle.write(','.join(reads))


if __name__ == "__main__":
    dump_read_ids(sys.argv[1])
