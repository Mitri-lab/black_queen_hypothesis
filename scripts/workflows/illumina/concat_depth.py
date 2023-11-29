import pandas as pd
from samples import Samples
from os.path import join
import sys

s = Samples()

contigs = {'ct':{'tig00000001_polypolish':'ct_0',
                 'tig00000002_polypolish':'ct_1'},
            'at':{'tig00000001_polypolish':'at_0',
                  'tig00000002_polypolish':'at_1',
                  'tig00000003_polypolish':'at_2',
                  'tig00000004_polypolish':'at_3',
                  'tig00000005_polypolish':'at_4'},
            'ms':{'tig00000001_polypolish':'ms_0'},
            'oa':{'tig00000002_polypolish':'oa_0',
                 'tig00000003_polypolish':'oa_1',
                 'tig00000004_polypolish':'oa_2',
                 'tig00000005_polypolish':'oa_3'}
                 }


def concat(df,strain):
      out = pd.DataFrame(columns=['chromosome', 'position', 'length'])

      i = -1
      prev_pos = 0
      prev_contig = None
      for contig, pos in zip(df['chromosome'], df['position']):
            if (prev_contig == contig) & (pos - 1 == prev_pos):
                  out.at[i, 'length'] += 1
            else:
                  i += 1
                  out.loc[i] = [contig, pos, 1]
            prev_pos = pos
            prev_contig = contig



      out['position'] = out['position'] - 2

      for i,row in out.iterrows():
            out.at[i,'chromosome'] = contigs[strain][row['chromosome']]
      return out


def caller(strain,sample):
      f = join(s.work,sample,strain,'depth_Q_0.tsv')
      df = pd.read_csv(f,sep='\t')
      df.columns = ['chromosome', 'position', 'coverage']
      df = df[df['coverage'] == 0]
      out = concat(df,strain)
      out.to_csv(join(s.work,sample,strain,'depth_Q_0.concat.csv'),index=False)
      return out
      

      
if __name__ == "__main__":
      strain = sys.argv[1]
      sample = sys.argv[2]
      out = caller(strain,sample)