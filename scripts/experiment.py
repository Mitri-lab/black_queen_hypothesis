from samples import Samples
from os.path import join
from os.path import exists
import pandas as pd

s = Samples()
class Experiment():
    def __init__(self):
        self.timepoints = ['T11','T22','T33','T44']
        self.treatments = [1,2,3,4]
        self.cosms = [1,2,3,4,5]
        self.strains = s.strains.keys()

        self.strain_series = {
                ('Agrobacterium tumefaciens',1):list(),
                ('Agrobacterium tumefaciens',3):list(),
                ('Agrobacterium tumefaciens',4):list(),
                ('Comamonas testosteroni',2):list(),
                ('Comamonas testosteroni',3):list(),
                ('Comamonas testosteroni',4):list(),
                ('Microbacterium saperdae',3):list(),
                ('Microbacterium saperdae',4):list(),
                ('Ochrobactrum anthropi',4):list()
                }

        for strain,samples in s.strains.items():
            for sample in samples:
                if sample['platform'] == 'illumina':
                    key = (strain,sample['treatment'])
                    self.strain_series[key].append(sample)