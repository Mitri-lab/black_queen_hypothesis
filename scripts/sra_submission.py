from samples import Samples
from os.path import join, exists
from os import sep, mkdir, removedirs, symlink, listdir
import pandas as pd
import subprocess


s = Samples()


ba = {i: i for i in s.strains.keys()}
ba["Ochrobactrum anthropi"] = "Brucella anthropi"
abbs = {"at": "At", "ct": "Ct", "oa": "Ba", "ms": "Ml"}

dir = join(
    sep,
    "work",
    "FAC",
    "FBM",
    "DMF",
    "smitri",
    "evomicrocomm",
    "genome_size",
    "data",
)

reps = {
    "T11.2.2": "T11.2.2.rep",
    "T11.2.4": "T11.2.4.rep",
    "T11.2.5": "T11.2.5.rep",
    "T33.2.1": "T33.2.1.rep",
    "T33.2.2": "T33.2.2.rep",
    "T33.2.3": "T33.2.3.rep",
    "T33.2.4": "T33.2.4.rep",
    "T33.2.5": "T33.2.5.rep",
    "T33.4.3": "T33.4.3.rep",
    "T44.2.2": "T44.2.2.rep",
    "T44.2.3": "T44.2.3.rep",
    "T44.2.4": "T44.2.4.rep",
    "T44.2.5": "T44.2.5.rep",
    "T44.3.3": "T44.3.3.rep",
    "T44.4.5": "T44.4.5.rep",
}
microcosms = {1: [1, 2], 2: [1, 2, 3, 4, 5], 3: [1, 2, 3, 4, 5], 4: [1, 2, 3, 4, 5]}
combinations = {
    1: "CAt",
    2: "CCt",
    3: "CAtCtMl",
    4: "CAtCtMlBa",
    0: "Ancestor",
}
timepoints = ["T11", "T22", "T33", "T44"]


def illumina():
    def create_read_links():
        for t in timepoints:
            for k, c in combinations.items():
                for m in microcosms[k]:
                    name = ".".join([t, c, "M" + str(m)])
                    target = join(dir, "SRA_submission", name)
                    if not exists(target):
                        mkdir(target)
                    folder = ".".join([t, str(k), str(m)])
                    if folder in reps.keys():
                        folder = reps[folder]
                    src_dir = join(dir, folder)
                    read1 = join(src_dir, "read1.fq.gz")
                    read2 = join(src_dir, "read2.fq.gz")
                    symlink(
                        read1, join(dir, "SRA_submission", target, name + "_1.fastq.gz")
                    )
                    symlink(
                        read2, join(dir, "SRA_submission", target, name + "_2.fastq.gz")
                    )

    def biosample():
        template = pd.read_csv(
            "Metagenome.environmental.1.0.tsv", sep="\t", comment="#"
        )
        out = pd.DataFrame(columns=template.columns)
        for t in timepoints:
            for k, c in combinations.items():
                for m in microcosms[k]:
                    name = ".".join([t, c, "M" + str(m)])
                    if k == 1:
                        species = "Agrobacterium tumefaciens"
                    elif k == 2:
                        species = "Comamonas testosteroni"
                    else:
                        species = "synthetic metagenome"
                    row = [
                        name,
                        "",
                        "PRJNA991498",
                        "synthetic metagenome",
                        "",
                        "industrial waste material [ENVO:00002267]",
                        "2014",
                        "United Kingdom:Melton Mowbray",
                        "52.76536031 N 0.87904024 W",
                        "https://pubmed.ncbi.nlm.nih.gov/15625673/",
                        "",
                        "",
                        "",
                        "",
                        "Week "
                        + t.split(".")[-1][1:]
                        + ", Microcosm "
                        + str(m)
                        + ", Combination "
                        + c[1:],
                        "",
                    ]
                    out.loc[len(out)] = row
        out.to_csv("BioSample.tsv", sep="\t", index=False)

    def metadata():
        template = pd.read_excel("SRA_metadata.xlsx", sheet_name="SRA_data")
        out = pd.DataFrame(columns=template.columns)
        i = 1
        for t in timepoints:
            for k, c in combinations.items():
                for m in microcosms[k]:
                    if (k == 1) | (k == 2):
                        library_source = "GENOMIC"
                    else:
                        library_source = "METAGENOMIC"
                    name = ".".join([t, c, "M" + str(m)])
                    name = ".".join([t, c, "M" + str(m)])
                    target = join(dir, "SRA_submission", name)
                    read1 = join(dir, "SRA_submission", target, name + "_1.fastq.gz")
                    read2 = join(dir, "SRA_submission", target, name + "_2.fastq.gz")
                    row = [
                        name,
                        i,
                        "Week "
                        + t.split(".")[-1][1:]
                        + ", Microcosm "
                        + str(m)
                        + ", Combination "
                        + c[1:],
                        "WGS",
                        library_source,
                        "RANDOM",
                        "paired",
                        "ILLUMINA",
                        "Illumina MiSeq",
                        "Populations of a microcosm were defrosted and regrown in metal working fluid. DNA was isolated using a Phenol-Chlorophorm-Isoamylalcohol solution.",
                        "fastq",
                        name + "_1.fastq.gz",
                        name + "_2.fastq.gz",
                        "",
                        "",
                        "",
                        "",
                    ]
                    i += 1
                    out.loc[len(out)] = row
        out.to_excel("metadata.xlsx", index=False)

    def ftp():
        for t in timepoints:
            for k, c in combinations.items():
                for m in microcosms[k]:
                    name = ".".join([t, c, "M" + str(m)])
                    target = join(dir, "SRA_submission", name)
                    read1 = join(dir, "SRA_submission", target, name + "_1.fastq.gz")
                    read2 = join(dir, "SRA_submission", target, name + "_2.fastq.gz")
                    cmd = ["sbatch", "transfer.sh", read1]
                    subprocess.call(" ".join(cmd), shell=True)
                    cmd = ["sbatch", "transfer.sh", read2]
                    subprocess.call(" ".join(cmd), shell=True)



def pacbio():
    def symlink():
        for strain, samples in s.strains.items():
            for sample in samples:
                if sample["platform"] == "pacbio":
                    src = join(sample["dir_name"], "reads.fastq.gz")
                    isolate = sample["name"][-1]
                    folder = ".".join(
                        [
                            sample["timepoint"],
                            abbs[s.abbreviations[strain]],
                            combinations[sample["treatment"]],
                            "M" + str(sample["cosm"]),
                            isolate,
                        ]
                    )
                    if not exists(join(dir, "SRA_submission", "pacbio", folder)):
                        mkdir(join(dir, "SRA_submission", "pacbio", folder))
                    trgt = join(
                        dir, "SRA_submission", "pacbio", folder, folder + ".fastq.gz"
                    )
                    symlink(src, trgt)

    def biosample():
        template = pd.read_csv("Microbe.1.0.tsv", sep="\t", comment="#")
        out = pd.DataFrame(columns=template.columns)
        for strain, samples in s.strains.items():
            for sample in samples:
                if sample["platform"] == "pacbio":
                    isolate = sample["name"][-1]
                    if strain == s.abbreviations["oa"]:
                        species = ba[strain]
                    elif strain == s.abbreviations["ms"]:
                        species = "Microbacterium liquefaciens"
                    else:
                        species = strain
                    name = ".".join(
                        [
                            sample["timepoint"],
                            abbs[s.abbreviations[strain]],
                            combinations[sample["treatment"]],
                            "M" + str(sample["cosm"]),
                            isolate,
                        ]
                    )
                    row = [
                        name,
                        "",
                        "PRJNA991498",
                        species,
                        "",
                        isolate,
                        "",
                        "Week "
                        + sample["timepoint"][-2:]
                        + ", Microcosm "
                        + str(sample["cosm"])
                        + ", Combination "
                        + combinations[sample["treatment"]]
                        + ", Isolate "
                        + isolate,
                        "2014",
                        "Switzerland: Lausanne",
                        "cell culture",
                        "",
                        "",
                        "",
                        "",
                        "",
                        "",
                        "",
                        "",
                        "",
                        "",
                        "",
                        "",
                        "Week " + sample["timepoint"][-2:] + ", serial transfers",
                        "",
                        "",
                        "",
                        "",
                        "",
                        "Week "
                        + sample["timepoint"][-2:]
                        + ", Microcosm "
                        + str(sample["cosm"])
                        + ", Combination "
                        + combinations[sample["treatment"]]
                        + ", Isolate "
                        + isolate,
                    ]
                    out.loc[len(out)] = row
        out.to_csv("BioSample.tsv", sep="\t", index=False)

    def metadata():
        template = pd.read_excel("SRA_metadata.xlsx", sheet_name="SRA_data")
        out = pd.DataFrame(columns=template.columns)
        i = 1
        out = pd.DataFrame(columns=template.columns)
        for strain, samples in s.strains.items():
            for sample in samples:
                if sample["platform"] == "pacbio":
                    isolate = sample["name"][-1]
                    name = ".".join(
                        [
                            sample["timepoint"],
                            abbs[s.abbreviations[strain]],
                            combinations[sample["treatment"]],
                            "M" + str(sample["cosm"]),
                            isolate,
                        ]
                    )
                    if strain == s.abbreviations["oa"]:
                        species = ba[strain]
                    elif strain == s.abbreviations["ms"]:
                        species = "Microbacterium liquefaciens"
                    else:
                        species = strain
                    target = join(dir, "SRA_submission", "pacbio", name)
                    read = join(
                        dir,
                        "SRA_submission",
                        "pacbio",
                        target,
                        name,
                        name + ".fastq.gz",
                    )
                    row = [
                        name,
                        i,
                        "Long-read sequencing of "
                        + species
                        + ", Isolate "
                        + isolate
                        + " of Microcosm "
                        + str(sample["cosm"])
                        + " of condition "
                        + str(combinations[sample["treatment"]]),
                        "WGS",
                        "GENOMIC",
                        "RANDOM",
                        "single",
                        "PACBIO_SMRT",
                        "Sequel II",
                        "Individual isolates from transfer 44 were sequenced.",
                        "fastq",
                        name + ".fastq.gz",
                        "",
                        "",
                        "",
                        "",
                        "",
                    ]
                    i += 1
                    out.loc[len(out)] = row
        out.to_excel("metadata.xlsx", index=False)

    def ftp():
        for strain, samples in s.strains.items():
            for sample in samples:
                if sample["platform"] == "pacbio":
                    isolate = sample["name"][-1]
                    name = ".".join(
                        [
                            sample["timepoint"],
                            abbs[s.abbreviations[strain]],
                            combinations[sample["treatment"]],
                            "M" + str(sample["cosm"]),
                            isolate,
                        ]
                    )
                    read = join(
                        dir, "SRA_submission", "pacbio", name, name + ".fastq.gz"
                    )
                    cmd = ["sbatch", "transfer.sh", read]
                    subprocess.call(" ".join(cmd), shell=True)


def symlink_rna():
    for strain, samples in s.strains.items():
        for sample in samples:
            if sample["platform"] == "RNA":
                src = join(sample["dir_name"])
                isolate = sample["name"].split(".")[-1]
                if sample["treatment"] == 0:
                    folder = ".".join(
                        [
                            "RNA",
                            "Ancestor",
                            abbs[s.abbreviations[strain]],
                            isolate,
                        ]
                    )
                else:
                    folder = ".".join(
                        [
                            "RNA",
                            sample["timepoint"],
                            abbs[s.abbreviations[strain]],
                            combinations[sample["treatment"]],
                            "M" + str(sample["cosm"]),
                            isolate,
                        ]
                    )
                if not exists(join(dir, "SRA_submission", "rna_seq", folder)):
                    mkdir(join(dir, "SRA_submission", "rna_seq", folder))
                trgt = join(
                    dir, "SRA_submission", "rna_seq", folder, folder + ".fastq.gz"
                )
                if not exists(trgt):
                    symlink(src, trgt)

def rna_seq():
    def biosample():
        template = pd.read_csv("Microbe.1.0.tsv", sep="\t", comment="#")
        out = pd.DataFrame(columns=template.columns)
        for strain, samples in s.strains.items():
            for sample in samples:
                if sample["platform"] == "RNA":
                    isolate = sample["name"].split(".")[-1]
                    if strain == s.abbreviations["oa"]:
                        species = ba[strain]
                    elif strain == s.abbreviations["ms"]:
                        species = "Microbacterium liquefaciens"
                    else:
                        species = strain
                    if sample["treatment"] == 0:
                        name = ".".join(
                            [
                                "RNA",
                                "Ancestor",
                                abbs[s.abbreviations[strain]],
                                isolate,
                            ]
                        )
                        week = "Ancestor"
                        isolation_source = "Ancestor " + ", Isolate " + isolate
                        description = "Ancestor"
                    else:
                        name = ".".join(
                            [
                                "RNA",
                                sample["timepoint"],
                                abbs[s.abbreviations[strain]],
                                combinations[sample["treatment"]],
                                "M" + str(sample["cosm"]),
                                isolate,
                            ]
                        )
                        week = "Week " + sample["timepoint"][-2:] + ", serial transfers"
                        isolation_source = (
                            "Week "
                            + sample["timepoint"][-2:]
                            + ", Microcosm "
                            + str(sample["cosm"])
                            + ", Combination "
                            + combinations[sample["treatment"]]
                            + ", Isolate "
                            + isolate
                        )
                        description = (
                            "Week "
                            + sample["timepoint"][-2:]
                            + ", Microcosm "
                            + str(sample["cosm"])
                            + ", Combination "
                            + combinations[sample["treatment"]]
                            + ", Isolate "
                            + isolate
                        )

                    row = [
                        name,
                        "",
                        "PRJNA991498",
                        species,
                        "",
                        isolate,
                        "",
                        isolation_source,
                        "2014",
                        "Switzerland: Lausanne",
                        "cell culture",
                        "",
                        "",
                        "",
                        "",
                        "",
                        "",
                        "",
                        "",
                        "",
                        "",
                        "",
                        "",
                        week,
                        "",
                        "",
                        "",
                        "",
                        "",
                        description,
                    ]
                    out.loc[len(out)] = row
        out.to_csv("BioSample.tsv", sep="\t", index=False)


    def metadata():
        template = pd.read_excel("SRA_metadata.xlsx", sheet_name="SRA_data")
        out = pd.DataFrame(columns=template.columns)
        i = 1
        out = pd.DataFrame(columns=template.columns)
        for strain, samples in s.strains.items():
            for sample in samples:
                if sample["platform"] == "RNA":
                    isolate = sample["name"].split(".")[-1]
                    if strain == s.abbreviations["oa"]:
                        species = ba[strain]
                    elif strain == s.abbreviations["ms"]:
                        species = "Microbacterium liquefaciens"
                    else:
                        species = strain
                    if sample["treatment"] == 0:
                        name = ".".join(
                            [
                                "RNA",
                                "Ancestor",
                                abbs[s.abbreviations[strain]],
                                isolate,
                            ]
                        )
                        title = "Short-read RNA sequencing of ancestral strain."
                        week = "Ancestor"
                        description = "The ancestor was sequenced from indiviual isolates."
                    else:
                        name = ".".join(
                            [
                                "RNA",
                                sample["timepoint"],
                                abbs[s.abbreviations[strain]],
                                combinations[sample["treatment"]],
                                "M" + str(sample["cosm"]),
                                isolate,
                            ]
                        )
                        title = "".join(
                            [
                                "Short-read RNA sequencing of "
                                + species
                                + ", Isolate "
                                + isolate
                                + " of Microcosm ",
                                str(sample["cosm"])
                                + " of condition "
                                + combinations[sample["treatment"]],
                            ]
                        )
                        week = ""
                        description = "Individual isolates from transfer 44 were sequenced."

                    target = join(dir, "SRA_submission", "pacbio", name)
                    read = join(
                        dir, "SRA_submission", "pacbio", target, name, name + ".fastq.gz"
                    )
                    row = [
                        name,
                        i,
                        title,
                        "RNA-Seq",
                        "TRANSCRIPTOMIC",
                        "cDNA",
                        "single",
                        "ILLUMINA",
                        "Illumina NovaSeq 6000",
                        description,
                        "fastq",
                        name + ".fastq.gz",
                        "",
                        "",
                        "",
                        "",
                        "",
                    ]
                    i += 1
                    out.loc[len(out)] = row
        out.to_excel("metadata.xlsx", index=False)

    def ftp():
        for strain, samples in s.strains.items():
            for sample in samples:
                if sample["platform"] == "RNA":
                    isolate = sample["name"].split(".")[-1]
                    if sample["treatment"] == 0:
                        name = ".".join(
                            [
                                "RNA",
                                "Ancestor",
                                abbs[s.abbreviations[strain]],
                                isolate,
                            ]
                        )

                    else:
                        name = ".".join(
                            [
                                "RNA",
                                sample["timepoint"],
                                abbs[s.abbreviations[strain]],
                                combinations[sample["treatment"]],
                                "M" + str(sample["cosm"]),
                                isolate,
                            ]
                        )
                    read = join(dir, "SRA_submission", "rna_seq", name, name + ".fastq.gz")
                    cmd = ["sbatch", "transfer.sh", read]
                    subprocess.call(" ".join(cmd), shell=True)