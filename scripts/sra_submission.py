from samples import Samples
from os.path import join, exists
from os import sep, mkdir, removedirs, symlink

s = Samples()

ba = {i: i for i in s.strains.keys()}
ba["Ochrobactrum anthropi"] = "Brucella anthropi"

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


def create_read_links():
    microcosms = {1: [1, 2], 2: [1, 2, 3, 4, 5], 3: [1, 2, 3, 4, 5], 4: [1, 2, 3, 4, 5]}
    combinations = {1: "CAt", 2: "CCt", 3: "CAtCtMl", 4: "CAtCtMlBa"}
    timepoints = ["T11", "T22", "T33", "T44"]

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
                symlink(read1,join(dir,'SRA_submission',target,name+'_1.fastq.gz'))
                symlink(read2,join(dir,'SRA_submission',target,name+'_2.fastq.gz'))

