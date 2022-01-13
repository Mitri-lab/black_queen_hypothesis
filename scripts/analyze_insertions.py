from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd
from os.path import join, exists
from samples import Samples
from plotting import Plotting
from hgt import Hgt

p = Plotting()
s = Samples()


def get_contigs(fasta):
    """Parses fastas and returns dictionary with contig name as
    key and sequence as value."""
    contigs = [contig for contig in SeqIO.parse(fasta, "fasta")]
    c = {contig.id: contig.seq for contig in contigs}
    return c


def filter_insertions(min_distance):
    """Filtering output from deleteion_detection.
    Filtering is necessary because many inserted sequences
    are located at start or end of contig which is probably
    a result of using different assemblers."""
    for strain, samples in s.strains.items():
        for sample in samples:
            if sample["platform"] == "pacbio":
                # Reading output of deletion detection
                df = pd.read_csv(
                    join(sample["dir_name"], "insertions.noalignments.tsv"),
                    sep="\t",
                ).drop_duplicates()
                # Adding insertion sequence for later processing
                df.insert(3, "sequence", None)
                fasta = join(sample["dir_name"], "assembly.fasta")
                contigs = get_contigs(fasta)
                # List which stors the index of items to delete because
                # they don't pass filter criteria
                to_pop = []
                # Iterating over every insertion
                for i, row in df.iterrows():
                    # Adding inserted sequence
                    df.at[i, "sequence"] = str(
                        contigs[row.chromosome][
                            row.position - 1 : row.position - 1 + row.length
                        ]
                    )
                    contig_length = len(contigs[row["chromosome"]])
                    # We keep insertion sequences which span the entire contig
                    if row["length"] == contig_length:
                        pass
                    # Filtering sequences at beginning of contig
                    elif row["position"] < min_distance:
                        to_pop.append(i)
                    # Filtering sequences at the end of contig
                    elif contig_length - row["position"] < min_distance:
                        to_pop.append(i)
                    else:
                        pass
                # Applying filter
                for element in to_pop:
                    df = df.drop(element)
                target = join(
                    sample["dir_name"], "insertions.noalignments.filtered.tsv"
                )
                # Dumping csv
                df.to_csv(target, sep="\t", index=False)


def track_insertions():
    """We iterate over every inserted nucleotide protein coding sequence
    and align it to every ancesteral strain present in the treatment.
    This allows us to see from where the insertions is coming.
    We were interested to see if the inserted sequences are potentially
    products of horizontal gene transfer.
    """
    for strain, samples in s.strains.items():
        for sample in samples:
            if sample["platform"] == "pacbio":
                hgt = Hgt(sample)
                f = join(sample["dir_name"], "insertions.noalignments.filtered.tsv")
                # Reading filtered inserted sequences
                df = pd.read_csv(f, sep="\t")
                df = df.dropna(subset=["nt_pos"])
                df.insert(len(df.columns), "origin", None)
                fasta = join(sample["dir_name"], "assembly.fasta")
                contigs = get_contigs(fasta)
                if len(df) > 0:
                    # Iterating over every inserted protein coding sequence
                    for i, row in df.iterrows():
                        id = row["chromosome"] + "." + str(row["position"])
                        # Position of protein sequence
                        gene_pos = [
                            int(element) for element in row["nt_pos"].split("-")
                        ]
                        # Getting sequence from assembly
                        seq = contigs[row["chromosome"]][gene_pos[0] : gene_pos[1]]
                        # Creating seq record
                        seq = SeqRecord(Seq(seq), id=id)
                        # This chunks the sequence with a window size of 500 and 250 steps
                        # Often this does nothing because sequence are often shorter than 500
                        chunks = hgt.chunker(seq, 500, 250)
                        # Aligns sequences to all present ancesteral genomes
                        hgt.mapper()
                        # Gets the contig name of the reference 
                        mapped_sequences = hgt.get_mapping_stats()
                        if len(mapped_sequences) > 0:
                            # Writing conting name of origin of sequence
                            df.at[i, "origin"] = " ".join(mapped_sequences)
                target = join(sample['dir_name'],'insertions.noalignments.filtered.tracked.tsv')
                # Dumping to csv         
                df.to_csv(target,sep='\t',index=False)
