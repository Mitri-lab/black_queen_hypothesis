from gc_bias import GC
import sys
from os.path import join

ref = sys.argv[1]
bam = sys.argv[2]
out = sys.argv[3]

gc = GC(ref, bam)
coverage = sum(list(gc.depth.values())) / len(gc.depth.values())
with open(join(out, 'coverage.txt'), 'w') as handle:
    handle.write(str(coverage))
