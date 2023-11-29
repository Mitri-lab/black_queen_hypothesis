from samples import Samples
from glob import glob
from os.path import join,split
import pandas as pd
import plotly.express as px

out = pd.DataFrame(columns=['strain','treatment','count'])

out 
s = Samples()
fs = glob(join(s.work,'rna_seq','ras_flow','MWF','genome','dea','countGroup','Ct*norm*'))
for f in fs:
    df = pd.read_csv(f,sep='\t',skiprows=1)
    strain,treatment = split(f)[-1][:2], split(f)[-1][2]
    i = 0
    for j,row in df.iterrows():
        if sum(list(row)[1:]) <= 2:
            i += 1
    out.loc[len(out)] = [strain,treatment,i]
    

