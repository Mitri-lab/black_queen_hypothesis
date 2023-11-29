from pysam import AlignmentFile
from samples import Samples
from os.path import join
import pandas as pd
import pysam
from subprocess import call, DEVNULL, STDOUT
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
s = Samples()

"""This is a script that we used to analyze a gene that we were interested.
Crp only shows up mutated in At evolved under condition 4.
It is suggested that could be a global transcription regulator.
Unfortunatly, we didn't get the RNA extraction to work for At.
"""

oa_seq = '/users/eulrich/work/genome_size/data/at/oa_sequence.fasta'
work = '/work/FAC/FBM/DMF/smitri/evomicrocomm/genome_size/data'

def get_crp_sequences():
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


def map_transposon():
    for strain, samples in s.strains.items():
        for sample in samples:
            if sample['platform'] == 'pacbio':
                reference = join(sample['dir_name'], 'assembly.fasta')
                sam = join(sample['dir_name'], 'oa_mapped.sam')
                cmd = [
                    "minimap2",
                    "-ax",
                    "asm5",
                    reference,
                    oa_seq,
                    ">",
                    sam,
                ]
                # Calling minimap and surpressing stdout
                call(" ".join(cmd), shell=True, stdout=DEVNULL,
                     stderr=STDOUT)


def print_stats():
    for strain, samples in s.strains.items():
        for sample in samples:
            if sample['platform'] == 'pacbio':
                a = pysam.AlignmentFile(join(sample['dir_name'], 'oa_mapped.sam'))
                counter = 0
                for read in a:
                    if not read.is_unmapped:
                        counter += 1
                print(sample['name'], counter)



def get_crp_sequences():
    crps = {'at': [('AGTU001.0001.c01', 'crp', 2351283, 2351738)],
        'ct': [('COTE001.0001.c01', 'crp_1', 2534712, 2535302),
               ('COTE001.0001.c01', 'crp_2', 4267275, 4268087)]}
    for strain, crps in crps.items():
        ref = join(work, strain, 'reference.fasta')
        contigs = {contig.id: contig for contig in SeqIO.parse(ref, 'fasta')}
        for crp in crps:
            contig = crp[0]
            name = crp[1]
            start = crp[2]
            end = crp[3]
            seq = contigs[contig][start:end]
            seq.id = name
            with open(join(work, strain, 'crps.fasta'), 'w') as handle:
                SeqIO.write(seq, handle, 'fasta')



def map_crps():
    for strain, samples in s.strains.items():
        for sample in samples:
            if (sample['platform'] == 'pacbio'):
                assembly = join(sample['dir_name'], 'assembly.fasta')
                crps = join(s.work, s.abbreviations[strain], 'crps.fasta')
                sam = join(sample['dir_name'], 'crps.sam')
                cmd = ['minimap2', '-ax', 'asm20', '--MD', assembly, crps, '>', sam]
                call(" ".join(cmd), shell=True, stdout=DEVNULL,
                    stderr=STDOUT)
