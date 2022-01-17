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

def get_contig_names():
    contig_names = {strain:None for strain in s.strains}
    reference_names = dict()
    for strain,reference in s.references.items():
        contig_names[strain] = [contig.id for contig in SeqIO.parse(reference,'fasta')]
        for contig_id in contig_names[strain]:
            reference_names[contig_id] = strain
    return contig_names,reference_names




sample = s.strains['Agrobacterium tumefaciens'][9]
contigs = get_contigs(join(sample['dir_name'],'assembly.fasta'))
origin = {contig:dict() for contig in contigs}
for name,contig in contigs.items():
    for position,base in enumerate(contig):
        origin[name][position] = sample['strain']
hgts = {key:dict() for key in contigs.keys()}
c_names,r_names = get_contig_names()
hgt = Hgt(sample)
window_size = 500
step = 100
hgts = []
for name,contig in contigs.items():
    record = SeqRecord(contig,id=name)
    seqs = hgt.chunker(record,window_size,step)
    hgt.mapper()
    mapped_seqs = hgt.get_mapping_stats()
    for name,hits in mapped_seqs.items():
        c_name = '.'.join(name.split('.')[:-1])
        pos = int(name.split('.')[-1])*step
        for ref,p,start,end,query_seq in hits:
            if ref not in c_names['Agrobacterium tumefaciens']:
                hgts.append((ref,p,query_seq))
                for j in range(pos+start,pos+end):
                    origin[c_name][j] = r_names[ref]