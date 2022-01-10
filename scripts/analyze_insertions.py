from Bio import SeqIO
import pandas as pd
from glob import glob
from os.path import join
from samples import Samples
from plotting import Plotting

p = Plotting()
s = Samples()


def get_contigs(fasta):
    """Parses fastas and returns dictionary with contig name as
    key and sequence as value."""
    contigs = [contig for contig in SeqIO.parse(fasta, "fasta")]
    c = {contig.id: contig.seq for contig in contigs}
    return c


def parse_insertions():
    """Parses all identified Pacbio sequences and stors them with
    inserted sequence in a dataframe."""
    dfs = []
    for strain, samples in s.strains.items():
        for sample in samples:
            if sample["platform"] == "pacbio":
                df = pd.read_csv(
                    join(sample["dir_name"], "mutant_to_parent.noalignments.tsv"),
                    sep="\t",
                    usecols=["chromosome", "position", "length"],
                ).drop_duplicates()
                df.insert(0,'strain',sample['strain'])
                df.insert(1, "sample_name", sample["name"])
                df.insert(2,'treatment',sample['treatment'])
                df.insert(6, "sequence", None)
                df.insert(7, "fasta", None)
                fasta = join(sample["dir_name"], "assembly.fasta")
                contigs = get_contigs(fasta)
                for i, row in df.iterrows():
                    df.at[i, "sequence"] = str(
                        contigs[row.chromosome][
                            row.position - 1 : row.position - 1 + row.length
                        ]
                    )
                    df.at[i, "fasta"] = contigs
                dfs.append(df)
    out = pd.concat(dfs)
    out.index = range(len(out))
    return out


def filter_insertions(insertions, min_distance):
    """To rule out assembly bias we look at insertion sequence
    which are not located at the end or the start of the assembly."""
    to_pop = []
    for i, row in insertions.iterrows():
        contig_length = len(row["fasta"][row["chromosome"]])
        if row["length"] == contig_length:
            pass
        elif row["position"] < min_distance:
            to_pop.append(i)
        elif contig_length - row["position"] < min_distance:
            to_pop.append(i)
        else:
            pass
    for element in to_pop:
        insertions = insertions.drop(element)
    return insertions

def plot_insertions(insertions):
    for strain in s.strains:
        df = insertions[insertions['strain'] == strain]
        samples = set(df['sample_name'])
        out = pd.DataFrame(0,columns=s.treatments[strain],index=samples)
        for sample,treatment,length in zip(df['sample_name'],df['treatment'],df['length']):
            out.at[sample,treatment] += length
        fig = p.subplot_treatments(strain,out)
        title = 'Filtered inserted bases in '+strain
        fig.update_layout(
            xaxis_title='samples',
            yaxis_title='inserted bp',
            title=title)
        fig.update_traces(showlegend=False)
        fig.write_image(join('..','plots','corrected_inserted_bases',title.replace(' ','_')+'.png'))



def get_sequences():
    work = "/users/eulrich/work/genome_size/data"
    df = pd.read_csv("sequences_of_interest.csv")
    df.insert(5, "sequence", None)
    for i, row in df.iterrows():
        c = row.chromosome
        p = row.position
        l = row.length
        s = row.sample_name
        fasta = join(work, s, "assembly.fasta")
        contigs = get_contigs(fasta)
        df.at[i, "sequence"] = str(contigs[c][p : p + l])
    return df
