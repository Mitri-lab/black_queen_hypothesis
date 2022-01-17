from os.path import join
from samples import Samples
from Bio import SeqIO
from subprocess import call
import subprocess
import pysam
import pandas as pd
from Bio.SeqRecord import SeqRecord


s = Samples()


class Hgt:
    """This class chunks a genetic sequence (can be a read or
    an assembly) with a sliding window algorithm. Chunks are
    then mapped to all ancestreal genomes from an experiment
    and the contig name of the reference where the chunks align
    are returned."""

    def __init__(self, sample):
        self.sample = sample
        self.seq_id = None

    def chunker(self, seq, window_size, step):
        """Creates chunks of a sequence. window_size defines
        chunk size and step the amount of basepair the windows
        is moved forward."""
        # List which stores all chunks
        seqs = []
        seqlen = len(seq)
        for counter, i in enumerate(range(0, seqlen, step)):
            # Returns ether entire sequence or window depending on sequence length
            j = seqlen if i + window_size > seqlen else i + window_size
            chunk = seq[i:j]
            self.seq_id = chunk.id
            # Add chunk id to sequence id
            chunk.id = chunk.id + "." + str(counter)
            seqs.append(chunk)
            if j == seqlen:
                break
        # Writes chunked sequence to hgt directory
        target = join(self.sample["dir_name"], "hgt",
                      "chunked_sequences.fasta")
        with open(target, "w") as handle:
            SeqIO.write(seqs, handle, "fasta")
        return seqs

    def mapper(self):
        """Maps chunked sequence to all ancesteral genomes
        present in experiment with minimap2.
        Minimap2 settings are set to accureate PacBio reads."""
        for strain in s.strains_per_treatment[self.sample["treatment"]]:
            reads = join(self.sample["dir_name"],
                         "hgt", "chunked_sequences.fasta")
            sam = join(self.sample["dir_name"], "hgt",
                       s.abbreviations[strain] + ".sam")
            cmd = [
                "minimap2",
                "-ax",
                "asm5",
                s.references[strain],
                reads,
                ">",
                sam,
            ]
            call(" ".join(cmd), shell=True, stdout=subprocess.DEVNULL,
                 stderr=subprocess.STDOUT)

    def get_mapping_stats(self):
        """Checks all mapped sequences and returns contig name
        of reference."""
        mapped_sequences = []
        mapped_sequences = dict()
        for strain in s.strains_per_treatment[self.sample["treatment"]]:
            sam = join(
                self.sample["dir_name"], "hgt", s.abbreviations[strain] + ".sam"
            )
            a = pysam.AlignmentFile(sam, "rb")
            reads = []
            # Iterating over all mapped primary reads
            for read in a:
                if (not read.is_unmapped) & (not read.is_secondary):
                    reads.append(read)
            # Appends contig name of reference
            for read in reads:
                if read.qname in mapped_sequences.keys():
                    mapped_sequences[read.qname].append(
                        (read.reference_name, read.pos,read.qstart,read.qend,SeqRecord(read.query_alignment_sequence,id=read.qname)))
                else:
                    mapped_sequences[read.qname] = []
                    mapped_sequences[read.qname].append(
                        (read.reference_name, read.pos,read.qstart,read.qend,SeqRecord(read.query_alignment_sequence,id=read.qname)))
        # Returns set of contig names of references
        return mapped_sequences
