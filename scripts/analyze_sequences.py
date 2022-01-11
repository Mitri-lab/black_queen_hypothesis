from os.path import join, exists
from os import mkdir, symlink
from samples import Samples
from Bio import SeqIO
from subprocess import call
import pysam

s = Samples()


def create_dirs():
    for ancestral, samples in s.strains.items():
        for sample in samples:
            if sample["platform"] == "pacbio":
                p = join(sample["dir_name"], "hgt")
                if not exists(p):
                    mkdir(p)
                for strain in s.strains_per_treatment[sample["treatment"]]:
                    reference = join(
                        sample["dir_name"], "hgt", s.abbreviations[strain] + ".fasta"
                    )
                    if not exists(reference):
                        symlink(
                            s.references[strain],
                            reference,
                        )


def chunker(seq, window_size, step):
    seqs = []
    seqlen = len(seq)
    for counter, i in enumerate(range(0, seqlen, step)):
        j = seqlen if i + window_size > seqlen else i + window_size
        chunk = seq[i:j]
        chunk.id = chunk.id + "." + str(counter)
        seqs.append(chunk)
        if j == seqlen:
            break
    return seqs


def mapper(sample):
    for strain in s.strains_per_treatment[sample["treatment"]]:
        reads = join(sample["dir_name"], "hgt", "chunked_reads.fasta")
        sam = join(sample["dir_name"], "hgt", s.abbreviations[strain] + ".sam")
        cmd = [
            "minimap2",
            "-ax",
            "asm5",
            s.references[strain],
            reads,
            ">",
            sam,
        ]
        call(" ".join(cmd), shell=True)


def get_cross_mapping(sample):
    cross_strains = [
        strain for strain in s.strains_per_treatment[sample['treatment']] if strain != sample["strain"]
    ]
    for cross_strain in cross_strains:
        sam = join(sample["dir_name"], "hgt", s.abbreviations[cross_strain] + ".sam")
        mapped_reads = []
        a = pysam.AlignmentFile(sam, "rb")
        for read in a:
            if (not read.is_unmapped) & (not read.is_secondary):
                mapped_reads.append(read)
                print(read.seq)


# sample = s.strains["Agrobacterium tumefaciens"][3]


def test():
    for strain, samples in s.strains.items():
        for sample in samples:
            if (sample["platform"] == "pacbio") & (sample["treatment"] not in [1,2]):
                print(sample["name"])
                fasta = join(sample["dir_name"], "corrected_reads.fasta")
                reads = [read for read in SeqIO.parse(fasta, "fasta")]
                chunks = []
                for read in reads:
                    chunks += chunker(read, 500, 500)
                target = join(sample["dir_name"], "hgt", "chunked_reads.fasta")
                with open(target, "w") as handle:
                    SeqIO.write(chunks, handle, "fasta")
                mapper(sample)
                get_cross_mapping(sample)
            
test()
