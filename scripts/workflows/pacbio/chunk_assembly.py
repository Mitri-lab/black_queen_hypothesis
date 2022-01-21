from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
from os.path import join


def chunk_assembly(fasta,out_dir):
    window_size = 1500
    step = 1500
    seqs = []
    contigs = [contig for contig in SeqIO.parse(fasta, 'fasta')]
    for contig in contigs:
        seqlen = len(contig)
        seq = contig
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
        target = join(out_dir, 
                      "chunked_assembly.fasta")
        with open(target, "w") as handle:
            SeqIO.write(seqs, handle, "fasta")

if __name__ == "__main__":
    chunk_assembly(sys.argv[1],sys.argv[2])