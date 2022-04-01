from pysam import AlignmentFile
from samples import Samples
from os.path import join
import pandas as pd
import pysam
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
s = Samples()
crps = []
names = ['At34.1', 'At44.1', 'At43.1', 'At45.1']
for strain, samples in s.strains.items():
    for sample in samples:
        if sample['name'] in names:
            print(sample['name'])
            f = join(sample['dir_name'], 'crps.sam')
            a = AlignmentFile(f)
            i = 0
            for crp in a:
                if not crp.is_unmapped:
                    name = sample['name'] + '.' + crp.qname + '.' + str(i)
                    print(name)
                    rec = SeqRecord(Seq(crp.get_reference_sequence()), id=name,name='crp')
                    if crp.is_reverse:
                        rec = rec.reverse_complement(id=name,name='crp')
                    crps.append(rec.translate(id=name,name='crp',description=name,to_stop=True))
                    i += 1
f = join(s.work, 'at', 'aa_crps.fasta')
with open(f, 'w') as handle:
    SeqIO.write(crps, handle, 'fasta')
