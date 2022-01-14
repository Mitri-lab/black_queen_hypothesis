from turtle import pos
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from samples import Samples
from os.path import join
from hgt import Hgt
import pandas as pd
from analyze_insertions import get_contigs

s = Samples()
def assembly():
    for strian,samples in s.strains.items():
        for sample in samples:
            if sample['platform'] == 'pacbio':
                hgt = Hgt(sample)
                fasta = join(sample['dir_name'],'assembly.fasta')
                reads = [read for read in SeqIO.parse(fasta,'fasta')]

                r_names = {key.id:None for key in reads}
                for read in reads:
                    seq = hgt.chunker(read,1500,500)
                    hgt.mapper()
                    r_names[read.id] = hgt.get_mapping_stats()

                contigs = []
                for names in r_names.values():
                    for name in names:
                        contigs.append(name)

                print(sample['name'],set(contigs))

def get_reference_names():
    r_names = {strain:None for strain in s.strains}
    for strain,reference in s.references.items():
        r_names[strain] = [contig.id for contig in SeqIO.parse(reference,'fasta')]
    return r_names



def reads():
    sample = s.strains['Agrobacterium tumefaciens'][9]
    contigs = get_contigs(join(sample['dir_name'],'assembly.fasta'))
    hgts = {key:dict() for key in contigs.keys()}
    reference = get_contigs(s.references[sample['strain']])
    r_names = get_reference_names()
    hgt = Hgt(sample)
    window_size = 50000
    step = 50000
    for name,contig in contigs.items():
        record = SeqRecord(contig,id=name)
        seqs = hgt.chunker(record,window_size,step)
        for seq in seqs:
            c = '.'.join(seq.id.split('.')[:-1])
            p = seq.id.split('.')[-1]*step
            target = join(sample['dir_name'],'hgt','chunked_sequences.fasta')
            SeqIO.write(seq,target,'fasta')
            hgt.mapper()
            positions = hgt.get_mapping_stats()
            for n,p in positions:
                if n not in r_names[sample['strain']]:
                    print(seq,positions)

    