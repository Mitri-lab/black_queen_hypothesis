from os.path import join, exists
from os import mkdir, symlink
from samples import Samples
from Bio import SeqIO
from subprocess import call
import subprocess
import pysam
import pandas as pd
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

s = Samples()


class Hgt:
    def __init__(self, sample):
        self.sample = sample
        self.seq_id = None

    def chunker(self, seq, window_size, step):
        seqs = []
        seqlen = len(seq)
        for counter, i in enumerate(range(0, seqlen, step)):
            j = seqlen if i + window_size > seqlen else i + window_size
            chunk = seq[i:j]
            self.seq_id = chunk.id
            chunk.id = chunk.id + "." + str(counter*step)
            seqs.append(chunk)
            if j == seqlen:
                break
        target = join(self.sample["dir_name"], "hgt", "chunked_sequences.fasta")
        with open(target, "w") as handle:
            SeqIO.write(seqs, handle, "fasta")

    def mapper(self):
        for strain in s.strains_per_treatment[self.sample["treatment"]]:
            reads = join(self.sample["dir_name"], "hgt", "chunked_sequences.fasta")
            sam = join(self.sample["dir_name"], "hgt", s.abbreviations[strain] + ".sam")
            cmd = [
                "minimap2",
                "-ax",
                "asm5",
                s.references[strain],
                reads,
                ">",
                sam,
            ]
            call(" ".join(cmd), shell=True,stdout=subprocess.DEVNULL,
    stderr=subprocess.STDOUT)

    def get_mapping_stats(self):
        mapped_sequences = []
        for strain in s.strains_per_treatment[self.sample["treatment"]]:
            sam = join(
                self.sample["dir_name"], "hgt", s.abbreviations[strain] + ".sam"
            )
            a = pysam.AlignmentFile(sam, "rb")
            reads = []
            for read in a:
                if (not read.is_unmapped) & (not read.is_secondary):
                    reads.append(read)
            for read in reads:
                mapped_sequences.append(read.reference_name)
        return set(mapped_sequences)
                 
                    


"""sample = s.strains["Agrobacterium tumefaciens"][9]
hgt = Hgt(sample)
f = join(
    sample["dir_name"], "mutant_to_parent.noalignments.filtered.tsv"
)
df = pd.read_csv(f, sep="\t")
if len(df) > 0:
    for i, row in df.iterrows():
        id = row["chromosome"] + "." + str(row["position"])
        seq = SeqRecord(Seq(row["sequence"].upper()), id=id, name=id)
        chunks = hgt.chunker(seq, 1000, 500)
        hgt.mapper()
        m = hgt.get_mapping_stats()
        break"""

