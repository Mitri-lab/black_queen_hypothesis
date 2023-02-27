import pysam
import numpy as np


f = '/users/eulrich/work/genome_size/data/rna_seq/ras_flow/tmp/MWF/genome/bamFileSort/10-1_At.3.3.1.1.sort.bam'

a = pysam.AlignmentFile(f)

coverage = {c: [] for c in a.references}
for c in a.references:
    counts = np.array(a.count_coverage(c))
    pos = range(a.get_reference_length(c))
    for p in pos:
        coverage[c].append(sum(counts[:, p]))
    break
