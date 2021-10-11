import re
import pandas as pd

class Samples():
    """This is a class used for parsing. It combines samples
    of Illumina and Pacbio"""
    def __init__(self):
        self.sample_sheet = './sample_sheet.csv'
        self.df = pd.read_csv(self.sample_sheet,dtype=str)
        self.strains = dict()

        #Adding reference fastas to class
        self.references = dict()
        self.references['Comamonas testosteroni'] = '/users/eulrich/evomicrocomm/references/ct/ct.fasta'
        self.references['Agrobacterium tumefaciens'] = '/users/eulrich/evomicrocomm/references/at/at.fasta'
        self.references['Microbacterium saperdae'] = '/users/eulrich/evomicrocomm/references/ms/ms.fasta'
        self.references['Ochrobactrum anthropi'] = '/users/eulrich/evomicrocomm/references/oa/oa.fasta'
      
        """For every strain we create a  list of dictionaries with
        the information about directory, sample name, treatment and sequencing platform"""
        strains = set(self.df['strain'])
        for strain in strains:
            self.strains[strain] = []
            df = self.df[self.df['strain'] == strain]
            for i,row in df.iterrows():
                meta = dict()
                meta['dir'] = row.dir_name
                meta['name'] = row.sample_name
                if row.platform == 'pacbio':
                    meta['treatment'] = int(row.sample_name[2])
                if row.platform == 'illumina':
                    meta['treatment'] = int(row.sample_name[4])
                meta['platform'] = row.platform
                self.strains[strain].append(meta)