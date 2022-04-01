from subprocess import call, DEVNULL, STDOUT
from samples import Samples
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from os.path import join
from pysam import AlignmentFile as af
s = Samples()

oa_seq = '/users/eulrich/work/genome_size/data/at/oa_sequence.fasta'
work = '/work/FAC/FBM/DMF/smitri/evomicrocomm/genome_size/data'
"""dirty code dump before clean implementation"""


"""
oa transposon stats
"""
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
                a = af(join(sample['dir_name'], 'oa_mapped.sam'))
                counter = 0
                for read in a:
                    if not read.is_unmapped:
                        counter += 1
                print(sample['name'], counter)



crps = {'at': [('AGTU001.0001.c01', 'crp', 2351283, 2351738)],
        'ct': [('COTE001.0001.c01', 'crp_1', 2534712, 2535302),
               ('COTE001.0001.c01', 'crp_2', 4267275, 4268087)]}

def get_crp_sequences():
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
