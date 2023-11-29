from samples import Samples
from os.path import join, exists,split
import pandas as pd
import plotly.express as px
from Bio import SeqIO
import vcfpy

"""This script hosts the code for the structural anlysis of the PacBio data."""

s = Samples()

strains = {s.abbreviations['at']: 'At',
           s.abbreviations['ct']: 'Ct',
           s.abbreviations['ms']: 'Ml',
           s.abbreviations['oa']: 'Oa'}

colors = {'ct': ['#7570B3', '#E6AB02', '#D95F02'],
          'at': ['#1B9E77', '#E6AB02', '#D95F02'],
          'all': ['#1B9E77', '#7570B3', '#E6AB02', '#D95F02'],
          'cosm': ['#4b2991', '#a431a0', '#ea4f88', '#f89178', '#edd9a3']
          }
colors_t = {'1': '#1B9E77',
                '2': '#7570B3',
                '3': '#E6AB02',
                '4': '#D95F02'}

# Execute parsing on cluster

w = 400
h = 250

def dump_pacbio_zero_cv():
    df = pd.DataFrame(columns=['strain', 'treatment', 'cosm', 'name','contig','start','end','length'])
    for strain, samples in s.strains.items():
        for sample in samples:
            if sample['platform'] == 'pacbio':
                f = join(sample['dir_name'],'depth_Q_0.concat.csv')
                if exists(f):
                    dels = pd.read_csv(f)
                    for i,row in dels.iterrows():
                        df.loc[len(df)] = [strain,sample['treatment'],sample['cosm'],sample['name'],row['chromosome'], \
                                        int(row['position']),int(row['position']) + int(row['length']),int(row['length'])]
    df.to_csv(join('..','variants','depth_Q_0.concat.csv'))


def parse_pacbio_deletions():
    df = pd.DataFrame(columns=['strain', 'treatment', 'cosm', 'deletions'])
    ins = pd.DataFrame(columns=['strain', 'treatment', 'cosm', 'insertions']) 
    for strain, samples in s.strains.items():
        for sample in samples:
            if sample['platform'] == 'pacbio':
                f = join(sample['dir_name'],'sniffles.var.vcf')
                if exists(f):
                    l = 0
                    k = 0
                    reader = vcfpy.Reader.from_path(f)
                    for record in reader:
                        for i, r in enumerate(record.ALT):
                            if record.INFO['SVTYPE'] == 'DEL':
                                l += record.INFO['SVLEN'] * (-1)
                            if record.INFO['SVTYPE'] == 'INS':
                                k += record.INFO['SVLEN'] 
                df.loc[len(df)] = [strain,sample['treatment'],sample['cosm'],l]
                ins.loc[len(ins)] = [strain,sample['treatment'],sample['cosm'],k]
    df.to_csv(join('..','variants','pacbio_inread_deletions.csv'),index=False)
    ins.to_csv(join('..','variants','pacbio_inread_insertions.csv'),index=False)

def parse_pac_no_cov():
    df = pd.DataFrame(columns=['strain', 'treatment', 'cosm', 'deletions'])
    i = 0
    species = [s.abbreviations[sp] for sp in ['at', 'ct','oa','ms']]
    for strain in species:
        for sample in s.strains[strain]:
            if sample['platform'] == 'pacbio':
                no_cov = pd.read_csv(join(sample['dir_name'],'depth_Q_0.concat.csv'))
                df.loc[len(df)] = [sample['strain'],sample['treatment'],sample['cosm'],sum(no_cov['length'])]
    df.to_csv(join('..', 'variants', 'pacbio_no_cov.csv'), index=False)

def dump_genes():
    df = pd.DataFrame(columns=['name','genes','treatment','strain','cosm'])
    #refs = [join(split(r),'bakta','assembly.contigs.polypolish.gbff') for r in s.references

    for strain,samples in s.strains.items():
        for sample in samples:
            if (sample['platform'] == 'pacbio'):
                f = join(sample['dir_name'],'bakta','assembly.contigs.racon.tsv')
                if exists(f):
                    tsv = pd.read_csv(f,sep='\t',skiprows=2)
                    tsv = tsv[tsv['Type'] == 'cds']
                    #tsv = tsv.dropna(subset=['Gene'])
                    df.loc[len(df)] = [sample['name'],len(tsv),sample['treatment'],strain,sample['cosm']]
                else:
                    pass

    df.to_csv(join('..','variants','genes.csv'),index=False)


def no_cov__illumina():
    # Broken, left if just in case it becomes relevant
    df = pd.DataFrame(columns=['strain', 'treatment',
                      'timepoint', 'cosm', 'deletions'])
    i = 0
    for species,samples in s.strains.items():
        for sample in samples:
            if sample['platform'] == 'illumina':
                f_0 = join(sample['dir_name'], 'depth_Q_0.tsv')
                f_60 = join(sample['dir_name'], 'depth_Q_60.tsv')
                fs = [f_0, f_60]
                for f in fs:
                    if exists(f):
                        ds = pd.read_csv(f, sep='\t', names=[
                                        'contig', 'pos', 'depth'])
                        mask = []
                        for c, d in zip(ds['contig'], ds['depth']):
                            if (c[:2] == s.abbreviations[strain]) & (d == 0):
                                mask.append(True)
                            else:
                                mask.append(False)
                        ds = ds[mask]
                        df.loc[i] = [strains[species], sample['treatment'],
                                    sample['timepoint'], sample['cosm'], len(ds)]
                        i += 1
    f = join('..', 'variants', 'deletions_illumina.csv')
    f_marc = join('..', 'variants', 'deletions_illumina_marc.csv')
    df.to_csv(f_marc)


def parse_zero_cov_in_read():
    df = pd.read_csv(join('..','variants','variants_comp_mapping.csv'))
    df = df[df['type'] == 'del']
    out = pd.DataFrame(columns=['strain','cosm','timepoint','treatment','deletions'])
    for strain,samples in s.strains.items():
        for sample in samples:
            if sample['platform'] == 'illumina':
                j = 0
                depth = pd.read_csv(join(sample['dir_name'],'depth_Q_0.concat.csv'))
                mask = (df['strain'] == strain) & (df['name'] == sample['name'])
                dels = df[mask]
                for i,row in depth.iterrows():
                    c, start, stop,length = \
                        row['chromosome'],row['position'],int(row['position']) + int(row['length']), row['length']
                    for i, d in dels[dels['chrom'] == c].iterrows():
                        if d['pos']  in [start,stop]:
                            j+= length
                            print(d)
                out.loc[len(out)] = strain, sample['cosm'], sample['timepoint'],sample['treatment'], j
                
                
def parse_assemblies():
    df = pd.DataFrame(columns=['strain', 'treatment', 'cosm', 'deletions','contigs'])
    i = 0
    species = [s.abbreviations[sp] for sp in ['at', 'ct','oa','ms']]
    for strain in species:
        for sample in s.strains[strain]:
            if sample['platform'] == 'pacbio':
                ds = 0  
                f = join(sample['dir_name'], 'assembly.contigs.racon.fasta')
                if exists(f):
                    contigs = [len(contig)
                            for contig in SeqIO.parse(f, 'fasta')]
                    l, n_contigs = sum(contigs), len(contigs)
                    
                    df.loc[i] = [strains[strain],
                                 sample['treatment'], sample['cosm'], l,n_contigs]
                    i += 1
    df.to_csv(join('..', 'variants', 'assembly_length.csv'), index=False)

# Local plotting


def font_size(fig, marker_size=3):
    """Style function for figures setting fot size and true black color."""
    for d in fig['data']:
        d['marker']['size'] = marker_size
        d['line']['width'] = 0.5
    # Font size
    j = 10
    fig.update_layout(font={'size': j, 'color': 'black'})
    for a in fig['layout']['annotations']:
        a['font']['size'] = j
        a['font']['color'] = 'black'
    fig['layout']['title']['font']['size'] = j
    fig['layout']['title']['font']['color'] = 'black'
    fig['layout']['legend']['title']['font']['size'] = j
    fig['layout']['legend']['title']['font']['color'] = 'black'
    fig.for_each_xaxis(lambda axis: axis.title.update(
        font=dict(size=j, color='black')))
    fig.for_each_yaxis(lambda axis: axis.title.update(
        font=dict(size=j, color='black')))
    fig.update_layout(
        margin=dict(l=20, r=20, t=20, b=20),
    )
    return fig

def plot_insertions_pacbio(abb):
    df = pd.read_csv(join('..','variants','pacbio_inread_insertions.csv'))
    mask = (df['strain'] == s.abbreviations[abb])
    df = df[mask]
    df = df.sort_values(by='treatment')
    df = df.replace(to_replace=0,value=0.1)
    fig = px.box(df, x='treatment', y='insertions', color='treatment', 
                    points='all', height=h, width=w/3,log_y=True)
    fig.update_traces(boxmean=True, quartilemethod="linear",
                        pointpos=0, jitter=1)
    fig.update_layout(showlegend=False, title='PacBio inread deletions', title_x=0.5)
    fig.update_xaxes(title='Condition',type='category')
    fig.update_yaxes(title='Deleted bp', rangemode='tozero')
    fig = font_size(fig,marker_size=5)
    fig.update_yaxes(title_standoff=0)
    fig.update_xaxes(title_standoff=0)
    for d in fig['data']:
        d['marker']['color'] = colors_t[d['name']]
        d['line']['color'] = colors_t[d['name']]
    fig.write_image(join('..', 'plots', 'plots', abb+'_pacbio_inread_insertions.svg'))

def plot_deletions_pacbio(abb):
    df = pd.read_csv(join('..','variants','pacbio_inread_deletions.csv'))
    mask = (df['strain'] == s.abbreviations[abb])
    df = df[mask]
    df = df.sort_values(by='treatment')
    #df = df.replace(to_replace=0,value=0.1)
    fig = px.box(df, x='treatment', y='deletions', color='treatment', 
                    points='all', height=h, width=w/3,log_y=False)
    fig.update_traces(boxmean=True, quartilemethod="linear",
                        pointpos=0, jitter=1)
    fig.update_layout(showlegend=False, title=abb + ' transfer 44', title_x=0.5)
    fig.update_xaxes(title='Condition',type='category')
    fig.update_yaxes(title='Within read deleted base pairs', rangemode='tozero')
    fig = font_size(fig,marker_size=5)
    fig.update_yaxes(title_standoff=0)
    fig.update_xaxes(title_standoff=0)
    for d in fig['data']:
        d['marker']['color'] = colors_t[d['name']]
        d['line']['color'] = colors_t[d['name']]
    fig.write_image(join('..','plots','plots',abb+'_pacbio_inread_deletions.svg'))

def plot_assembly_length(abb):
    genome_length = {'ct':6130738,
                    'at':5430046,
                    'oa':5249884,
                    'ms':3761082}

    """ancestors = ['ct','at']
    for a in ancestors:
        l = 0 
        f = join(s.work,'ancestors',a,'assembly.contigs.polypolish.fasta')
        contigs = [contig for contig in SeqIO.parse(f,'fasta')]
        for c in contigs:
            l += len(c)
            genome_length[a] = l"""


    df = pd.read_csv(join('..', 'variants', 'assembly_length.csv'))
    df = df.replace(to_replace='Ct',value='ct')
    df = df.replace(to_replace='At',value='at')
    df = df.replace(to_replace='Oa',value='oa')
    df = df.replace(to_replace='Ml',value='ms')
    mask = (df['strain'] == abb)
    df = df[mask]
    df = df.sort_values(by='treatment')
    #df = df.replace(to_replace=0,value=0.1)
    hover_data = ['treatment','cosm']
    fig = px.box(df, x='treatment', y='deletions', color='treatment', 
                    points='all', height=h, width=w/3,log_y=True,hover_data=hover_data)
    fig.add_hline(genome_length[abb],line_dash="dot",line_color='gray',opacity=0.75)
    fig.update_traces(boxmean=True, quartilemethod="linear",
                        pointpos=0, jitter=1)
    fig.update_layout(showlegend=False, title=abb + ' transfer 44', title_x=0.5)
    fig.update_xaxes(title='Condition',type='category')
    fig.update_yaxes(title='Assebmly length [base pairs]', rangemode='tozero')
    fig = font_size(fig,marker_size=5)
    fig.update_yaxes(title_standoff=0)
    fig.update_xaxes(title_standoff=0)
    for d in fig['data']:
        d['marker']['color'] = colors_t[d['name']]
        d['line']['color'] = colors_t[d['name']]
    fig.write_image(join('..','plots','plots',abb+'_assembly_length.svg'))
    return fig

def plot_contigs(abb):
    genome_length = {'ct':2,
                    'at':5,
                    'oa':4,
                    'ms':1}

    """ancestors = ['ct','at']
    for a in ancestors:
        l = 0 
        f = join(s.work,'ancestors',a,'assembly.contigs.polypolish.fasta')
        contigs = [contig for contig in SeqIO.parse(f,'fasta')]
        for c in contigs:
            l += len(c)
            genome_length[a] = l"""


    df = pd.read_csv(join('..', 'variants', 'assembly_length.csv'))
    df = df.replace(to_replace='Ct',value='ct')
    df = df.replace(to_replace='At',value='at')
    df = df.replace(to_replace='Oa',value='oa')
    df = df.replace(to_replace='Ml',value='ms')

    mask = (df['strain'] == abb)
    df = df[mask]
    df = df.sort_values(by='treatment')
    #df = df.replace(to_replace=0,value=0.1)
    hover_data = ['treatment','cosm']
    fig = px.box(df, x='treatment', y='contigs', color='treatment', 
                    points='all', height=h, width=w/3,log_y=True,hover_data=hover_data)
    fig.add_hline(genome_length[abb],line_dash="dot",line_color='gray',opacity=0.75)
    fig.update_traces(boxmean=True, quartilemethod="linear",
                        pointpos=0, jitter=1)
    fig.update_layout(showlegend=False, title=abb + ' transfer 44', title_x=0.5)
    fig.update_xaxes(title='Condition',type='category')
    fig.update_yaxes(title='N contigs', rangemode='tozero')
    fig = font_size(fig,marker_size=5)
    fig.update_yaxes(title_standoff=0)
    fig.update_xaxes(title_standoff=0)
    for d in fig['data']:
        d['marker']['color'] = colors_t[d['name']]
        d['line']['color'] = colors_t[d['name']]
    fig.write_image(join('..','plots','plots',abb+'_contigs.svg'))


def pacbio_zero_coverage(abb):
    df = pd.read_csv(join('..', 'variants', 'pacbio_no_cov.csv'))
    mask = (df['strain'] == s.abbreviations[abb])
    df = df[mask]
    df = df.sort_values(by='treatment')
    #df = df.replace(to_replace=0,value=0.1)
    hover_data = ['treatment','cosm']
    fig = px.box(df, x='treatment', y='deletions', color='treatment', 
                    points='all', height=h, width=w/3,log_y=True,hover_data=hover_data)
    fig.update_traces(boxmean=True, quartilemethod="linear",
                        pointpos=0, jitter=1)
    fig.update_layout(showlegend=False, title=abb + ' transfer 44', title_x=0.5)
    fig.update_xaxes(title='Condition',type='category')
    fig.update_yaxes(title='Base pairs with 0 PacBio coverage', rangemode='tozero')
    fig = font_size(fig,marker_size=5)
    fig.update_yaxes(title_standoff=0)
    fig.update_xaxes(title_standoff=0)
    for d in fig['data']:
        d['marker']['color'] = colors_t[d['name']]
        d['line']['color'] = colors_t[d['name']]
    fig.write_image(join('..','plots','plots',abb+'_pacbio_zero_coverage.svg'))
    return fig

def number_of_genes(abb):
    df = pd.read_csv(join('..', 'variants', 'genes.csv'))
    mask = (df['strain'] == s.abbreviations[abb])
    df = df[mask]
    df = df.sort_values(by='treatment')
    #df = df.replace(to_replace=0,value=0.1)
    hover_data = ['treatment','cosm']
    fig = px.box(df, x='treatment', y='genes', color='treatment', 
                    points='all', height=h, width=w/3,log_y=True,hover_data=hover_data)
    fig.update_traces(boxmean=True, quartilemethod="linear",
                        pointpos=0, jitter=1)
    fig.update_layout(showlegend=False, title=abb + ' transfer 44', title_x=0.5)
    fig.update_xaxes(title='Condition',type='category')
    fig.update_yaxes(title='Number of coding sequences', rangemode='tozero')
    fig = font_size(fig,marker_size=5)
    fig.update_yaxes(title_standoff=0)
    fig.update_xaxes(title_standoff=0)
    for d in fig['data']:
        d['marker']['color'] = colors_t[d['name']]
        d['line']['color'] = colors_t[d['name']]
    fig.write_image(join('..','plots','plots',abb+'_genes.svg'))
