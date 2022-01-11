from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd
from os.path import join, exists
from samples import Samples
from plotting import Plotting
from analyze_sequences import Hgt

p = Plotting()
s = Samples()


def get_contigs(fasta):
    """Parses fastas and returns dictionary with contig name as
    key and sequence as value."""
    contigs = [contig for contig in SeqIO.parse(fasta, "fasta")]
    c = {contig.id: contig.seq for contig in contigs}
    return c


def fiilter_insertions(min_distance):
    dfs = []
    for strain, samples in s.strains.items():
        for sample in samples:
            if sample["platform"] == "pacbio":
                df = pd.read_csv(
                    join(sample["dir_name"], "mutant_to_parent.noalignments.tsv"),
                    sep="\t",
                    usecols=["chromosome", "position", "length"],
                ).drop_duplicates()
                df.insert(3, "sequence", None)
                fasta = join(sample["dir_name"], "assembly.fasta")
                contigs = get_contigs(fasta)
                to_pop = []
                for i, row in df.iterrows():
                    df.at[i, "sequence"] = str(
                        contigs[row.chromosome][
                            row.position - 1 : row.position - 1 + row.length
                        ]
                    )
                    contig_length = len(contigs[row["chromosome"]])
                    if row["length"] == contig_length:
                        pass
                    elif row["position"] < min_distance:
                        to_pop.append(i)
                    elif contig_length - row["position"] < min_distance:
                        to_pop.append(i)
                    else:
                        pass
                for element in to_pop:
                    df = df.drop(element)
                target = join(
                    sample["dir_name"], "mutant_to_parent.noalignments.filtered.tsv"
                )
                df.to_csv(target, sep="\t")


def analyze_insertions():
    for strain, samples in s.strains.items():
        for sample in samples:
            if sample["platform"] == "pacbio":
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
                        target = join(
                            sample["dir_name"], "hgt", "chunked_sequences.fasta"
                        )
                        with open(target, "w") as handle:
                            SeqIO.write(chunks, handle, "fasta")
                        # hgt.mapper()
                        if hgt.get_cross_mapping():
                            print(
                                sample["name"],
                                row["chromosome"],
                                row["position"],
                                row["length"],
                                row["sequence"],
                            )


def plot_filtered_insertions():
    """This function plots the sum of deleted and inserted base pairs."""
    i = {strain: None for strain in s.strains}
    for strain in s.strains:
        treatments = s.treatments[strain]
        inserted_bases = pd.DataFrame(
            columns=treatments,
            index=[
                sample["name"]
                for sample in s.strains[strain]
                if sample["platform"] == "pacbio"
            ],
        )
        for sample in s.strains[strain]:
            if sample["platform"] == "pacbio":
                insertions = join(
                    sample["dir_name"], "mutant_to_parent.noalignments.filtered.tsv"
                )
                if exists(insertions):
                    inserted_bases.at[sample["name"], sample["treatment"]] = sum(
                        pd.read_csv(
                            insertions,
                            sep="\t",
                            usecols=["chromosome", "position", "length"],
                        ).drop_duplicates()["length"]
                    )
        i[strain] = inserted_bases
        fig = p.subplot_treatments(strain, inserted_bases)
        title = "Filtered inserted bases in " + strain
        fig.update_layout(xaxis_title="samples", yaxis_title="inserted bp", title=title)
        fig.update_traces(showlegend=False)
        fig.write_image(
            join(
                "..",
                "plots",
                "corrected_inserted_bases",
                title.replace(" ", "_") + ".png",
            )
        )

    return i
