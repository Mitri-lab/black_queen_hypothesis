import pandas as pd
from samples import Samples
from os.path import join
import sys

s = Samples()



def concat(df):
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



      out['position'] = out['position'] - 1

      return out


def caller(sample):
      f = join(s.work,sample,'depth_Q_0.tsv')
      df = pd.read_csv(f,sep='\t')
      df.columns = ['chromosome', 'position', 'coverage']
      df = df[df['coverage'] == 0]
      out = concat(df)
      out.to_csv(join(s.work,sample,'depth_Q_0.concat.csv'),index=False)
      return out
      

      
if __name__ == "__main__":
      sample = sys.argv[1]
      out = caller(sample)