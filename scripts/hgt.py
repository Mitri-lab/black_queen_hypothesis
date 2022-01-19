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
        # Dictionary with sample info from Sample class
        self.sample = sample
        # Fasta of sample assembly
        fasta = join(self.sample['dir_name'], 'assembly.fasta')
        # Contigs in dict form
        self.contigs = self.get_contigs(fasta)
        # Step size for sliding windows algorithm
        self.step = None
        # Dictionary storing origins of sequence in assembly
        self.origins = {key: dict() for key in self.contigs.keys()}
        # Filtered origins
        self.filtered = {key: dict() for key in self.contigs.keys()}

    def get_contig_names(self):
        """Returns dictionary. contig_names with strains as keys
        and contig names as values. reference_names with contig
        names as keys and strain as value."""
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
        chunk size and step the amount of basepairs the window
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
        return seqs

    def chunk_assembly(self):
        """Chunks an assembly of multiple contigs into different 
        chunks using a sliding window algorightm (see chunker function)."""
        assembly_chunks = []
        for name, contig in self.contigs.items():
            record = SeqRecord(contig, id=name)
            # Creates chunks of every contig
            assembly_chunks += self.chunker(record, 500, 100)
        target = join(self.sample["dir_name"], "hgt",
                      "chunked_sequences.fasta")
        # Dumps chunks to fasta
        with open(target, "w") as handle:
            SeqIO.write(assembly_chunks, handle, "fasta")

    def mapper(self):
        """Maps chunked sequence to all ancesteral genomes
        experiment with minimap2.
        Minimap2 settings are set to accureate PacBio reads."""
        for strain in s.strains.keys():
            # Query sequences created with chunk_assembly()
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
            # Calling minimap and surpressing stdout
            call(" ".join(cmd), shell=True, stdout=subprocess.DEVNULL,
                 stderr=subprocess.STDOUT)

    def get_mapping_stats(self):
        """Checks all mapped sequences and returns contig name
        of reference."""
        c_names, r_names = self.get_contig_names()
        for strain in s.strains.keys():
            # Iterating over every sam file
            sam = join(
                self.sample["dir_name"], "hgt", s.abbreviations[strain] + ".sam"
            )
            a = pysam.AlignmentFile(sam, "rb")
            reads = []
            # Iterating over all reads
            # Read must be primary,mapped and have quality of 60
            for read in a:
                if (not read.is_unmapped) & (not read.is_secondary) & (read.mapq == 60):
                    reads.append(read)
            # Appends contig name of reference
            for read in reads:
                # Contig and position of query sequence form assembly are
                # stored in contig name
                # First part is contig name, second number is n step
                # step size is known and position therefore as well
                name = '.'.join(read.qname.split('.')[:-1])
                pos = int(read.qname.split('.')[-1])*self.step
                # Iteratign over aligned query sequence
                for j in range(pos+read.query_alignment_start, pos+read.query_alignment_end):
                    # Checking if position already in dictionary
                    if j not in self.origins[name].keys():
                        self.origins[name][j] = []
                        # Appending strain at given position
                        self.origins[name][j].append(
                            r_names[read.reference_name])

                    else:
                        # Appending strain at given position
                        self.origins[name][j].append(
                            r_names[read.reference_name])

    def dump_origins(self):
        """Filters identified origins and dumps to json. Only positions
        with either two different origins or an origin which is not anceteral
        are of interest"""
        for name, contig in self.origins.items():
            for pos, strains in contig.items():
                # Applying filter
                if (len(set(strains)) > 1) or (strains[0] != self.sample['strain']):
                    self.filtered[name][pos] = list(set(strains))

        # Dumping to json
        j_f = json.dumps(self.filtered, indent=4)
        with open(join(self.sample['dir_name'], 'origins.json'), 'w') as handle:
            handle.write(j_f)
