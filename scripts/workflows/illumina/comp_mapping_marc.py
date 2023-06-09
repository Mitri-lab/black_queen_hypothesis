from samples import Samples
from Bio import SeqIO
from os.path import join
from subprocess import call
import pysam
import sys

s = Samples()

"""Maps reads to merged references. Competitive mapping allowing
only primary mapping. Bam files are filtered for reads with mapped in 
proper pairs, before bam files are split per species."""

def get_sample(name):
    for strain, samples in s.strains.items():
        for sample in samples:
            if sample['name'] == name:
                return sample


def get_ref_species(t):
    species = []
    for specie, treatments in s.treatments.items():
        if t in treatments:
            species.append(specie)
    return species


def dump_reference(name):
    s_trgt = get_sample(name)
    t = s_trgt['treatment']
    species = get_ref_species(t)
    refs = [join(s.work,s.abbreviations[sp],'assembly.fasta') for sp in species]

    out = []
    for r, sp in zip(refs, species):
        contigs = [contig for contig in SeqIO.parse(r, 'fasta')]
        for j, c in enumerate(contigs):
            name = s.abbreviations[sp] + '_' + str(j)
            c.id, c.name, c.description = name, name, ''
            out.append(c)

    with open(join(s.work, s_trgt['name'], 'reference_marc.fasta'), 'w') as handle:
        SeqIO.write(out, handle, 'fasta')


def map_reads(name):
    s_trgt = get_sample(name)
    out_f = join(s.work, s_trgt['name'], 'mapped_reads_marc')
    ref = join(s.work, s_trgt['name'], 'reference_marc.fasta')
    r1 = join(s.work, name, 'read1.fq.gz')
    r2 = join(s.work, name, 'read2.fq.gz')
    minimap = ['minimap2', '--secondary=no', '-t', '8', '-R',
               "'@RG\\tID:snippy\\tSM:snippy'", '-ax', 'sr', ref, r1, r2, '>', out_f + '.sam']
    call(' '.join(minimap), shell=True)
    samtools = ['samtools', 'view', '-b', '-f', '3', '-q', '60', out_f +
                '.sam', '|', 'samtools', 'sort', '--threads', '8',  '-o', out_f + '.sorted.bam']
    call(' '.join(samtools), shell=True)
    index = ['samtools', 'index', out_f + '.sorted.bam']
    call(' '.join(index), shell=True)


def split_bam(name):
    b = join(s.work, name, 'mapped_reads_marc.sorted.bam')
    a = pysam.AlignmentFile(b)
    s_trgt = get_sample(name)
    species = [s.abbreviations[sp]
               for sp in get_ref_species(s_trgt['treatment'])]
    outs = {sp: pysam.AlignmentFile(
        join(s.work, name, sp, 'mapped_reads_marc.bam'), 'wb', template=a) for sp in species}
    for read in a:
        if (not read.is_secondary):
            outs[read.reference_name[:2]].write(read)


def terminate(name):
    f = join(s.work, name, 'done_marc.txt')
    with open(f, 'w') as handle:
        handle.write(f)


if __name__ == "__main__":
    name = sys.argv[1]
    dump_reference(name)
    map_reads(name)
    split_bam(name)
    terminate(name)
