from os.path import join
from samples import Samples
from Bio import SeqIO
from subprocess import call
import subprocess
import pysam
from Bio.SeqRecord import SeqRecord
import json

s = Samples()


class Hgt:
    """This class chunks a genetic sequence (can be a read or
    an assembly) with a sliding window algorithm. Chunks are
    then mapped to all ancestreal genomes from an experiment
    and the contig name of the reference where the chunks align
    are returned."""

    def __init__(self, sample):
        self.sample = sample
        fasta = join(self.sample['dir_name'], 'assembly.fasta')
        self.contigs = self.get_contigs(fasta)

    def get_contig_names(self):
        contig_names = {strain: None for strain in s.strains}
        reference_names = dict()
        for strain, reference in s.references.items():
            contig_names[strain] = [
                contig.id for contig in SeqIO.parse(reference, 'fasta')]
            for contig_id in contig_names[strain]:
                reference_names[contig_id] = strain
        return contig_names, reference_names

    def get_contigs(self, fasta):
        """Parses fastas and returns dictionary with contig name as
        key and sequence as value."""
        contigs = [contig for contig in SeqIO.parse(fasta, "fasta")]
        c = {contig.id: contig.seq for contig in contigs}
        return c

    def chunker(self, seq, window_size, step):
        """Creates chunks of a sequence. window_size defines
        chunk size and step the amount of basepair the windows
        is moved forward."""
        # List which stores all chunks
        seqs = []
        seqlen = len(seq)
        self.step = step
        for counter, i in enumerate(range(0, seqlen, step)):
            # Returns ether entire sequence or window depending on sequence length
            j = seqlen if i + window_size > seqlen else i + window_size
            chunk = seq[i:j]
            # Add chunk id to sequence id
            chunk.id = chunk.id + "." + str(counter)
            seqs.append(chunk)
            if j == seqlen:
                break
        # Writes chunked sequence to hgt directory
        return seqs

    def chunk_assembly(self):
        assembly_chunks = []
        for name, contig in self.contigs.items():
            record = SeqRecord(contig, id=name)
            assembly_chunks += self.chunker(record, 1000, 500)
        target = join(self.sample["dir_name"], "hgt",
                      "chunked_sequences.fasta")
        with open(target, "w") as handle:
            SeqIO.write(assembly_chunks, handle, "fasta")
        self.chunks = assembly_chunks

    def mapper(self):
        """Maps chunked sequence to all ancesteral genomes
        present in experiment with minimap2.
        Minimap2 settings are set to accureate PacBio reads."""
        for strain in s.strains.keys():
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
        mapped_sequences = {key: dict() for key in self.contigs.keys()}
        c_names, r_names = self.get_contig_names()
        for strain in s.strains.keys():
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
                name = '.'.join(read.qname.split('.')[:-1])
                pos = int(read.qname.split('.')[-1])*self.step
                for j in range(pos+read.query_alignment_start, pos+read.query_alignment_end):
                    if j not in mapped_sequences[name].keys():
                        mapped_sequences[name][j] = []
                        mapped_sequences[name][j].append(
                            r_names[read.reference_name])

                    else:
                        mapped_sequences[name][j].append(
                            r_names[read.reference_name])
        # Returns set of contig names of references
        self.mapped_sequences = mapped_sequences

    def dump_origins(self):
        r_names,c_names = self.get_contig_names()
        self.filtered = {key: dict() for key in self.contigs.keys()}
        for name, contig in self.mapped_sequences.items():
            for pos, strains in contig.items():
                if (len(set(strains)) > 1) or (strains[0] != self.sample['strain']):
                    self.filtered[name][pos] = list(set(strains))

        j_f = json.dumps(self.filtered, indent=4)
        with open(join(self.sample['dir_name'], 'origins.json'), 'w') as handle:
            handle.write(j_f)
