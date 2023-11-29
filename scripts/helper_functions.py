from samples import Samples
import pandas as pd
from os.path import join,exists
from Bio import SeqIO
import re



"""
################################################################################
Author: https://github.com/nahanoo
This is a collection of functions for illumina data that are not relevant for the paper
but worth keeping.
################################################################################
"""

s = Samples()


def plot_effects():
    """This plots the summed effects grouped per treatment per species
    """
    snps = {strain: None for strain in s.strains}
    for strain, samples in s.strains.items():
        effects = []
        for sample in samples:
            if (sample['platform'] == 'illumina') & (sample['timepoint'] == 'T44'):
                f = join(sample['dir_name'], 'snippy', 'snps.tab')
                df = pd.read_csv(f, sep='\t').drop_duplicates()
                df = df[df['EFFECT'].notna()]
                for effect in df['EFFECT']:
                    effects.append(effect.split(' ')[0])
        effects = list(set(effects))
        columns = s.treatments[strain]
        snp = pd.DataFrame(columns=columns, index=effects)
        for sample in samples:
            if (sample['platform'] == 'illumina') & (sample['timepoint'] == 'T44'):
                f = join(sample['dir_name'], 'snippy', 'snps.tab')
                df = pd.read_csv(f, sep='\t').drop_duplicates()
                df = df[df['EFFECT'].notna()]
                effects = []
                for effect in df['EFFECT']:
                    effects.append(effect.split(' ')[0])

                for effect in set(effects):
                    mask = []
                    for e in df['EFFECT']:
                        if re.search(effect, e, flags=re.IGNORECASE):
                            mask.append(True)
                        else:
                            mask.append(False)
                    if pd.isna(snp.at[effect, sample['treatment']]):
                        snp.at[effect, sample['treatment']] = len(df[mask])
                    else:
                        snp.at[effect, sample['treatment']] += len(df[mask])
        snps[strain] = snp
    return snps


def write_gc_content():
    """Small function to get gc content of references."""
    def get_gc_content(sequence):
        """This function returns gc content for string sequence"""
        return 100.0*len([base for base in sequence if base in "GC"])/len(sequence)
    """This is a little helper function to get the GC content of the wild-type genomes"""
    # Setting up df and iterating over all strains
    df = pd.DataFrame(columns=['gc content'])
    for strain in s.strains:
        sequence = str()
        reference = s.references[strain]
        contigs = [contig for contig in SeqIO.parse(reference, 'fasta')]
        for contig in contigs:
            # Adding sequence as string
            sequence += contig.seq
        # Getting gc content of sequence
        gc_content = get_gc_content(sequence)
        # Writing to df and to file
        df.at[strain, 'gc content'] = gc_content
    df.index.name = 'strain'
    fname = join('..', 'tables', 'gc_content', 'gc_content.csv')
    df.to_csv(fname)

"""Functions to analyze GC content of deletions"""
def get_gc_content(sequence):
    """This function returns gc content for string sequence"""
    return 100.0*len([base for base in sequence if base in "GC"])/len(sequence)

def get_gc_contents():
    """This function stores the GC contents of every 150-mer of the reference sequences."""
    window_size = 150
    gc_contents = dict()
    for strain in s.strains:
        gc_content = []
        #Parsing reference sequence by making use of Samples class
        reference = {contig.name:contig for contig in SeqIO.parse(s.references[strain],'fasta')}
        #Iterating over every contig
        for chromosome,record in reference.items():
            #Iterating over every position of the contig and caclulating gc content for window.
            for position in range(len(record)-window_size):
                gc_content.append(get_gc_content(record[position:position+window_size]))
        gc_contents[strain] = gc_content
    return gc_contents

def gc_histogram(fig,gc_contents,strain):
    """Takes a plotly subplot as input and plots gc content as histogram."""
    fig.append_trace(go.Histogram(x=gc_contents[strain],name='gc content',nbinsx=20),1,1)
    fig.update_layout(
        xaxis_title = 'GC content in %',
        yaxis_title = 'counts'
    )
    return fig

def deletion_histogram(fig,strain):
    """Parses all deletions identified in Illumina. Calculates the GC
    content of the location of the deletion. Since many deletions are very small
    50 basepairs are added to both sides."""
    #Data path
    gc_contents = []
    for sample in s.strains[strain]:
        #Not all samples have deletions thus the os.path.exists check
        f = os.path.join(sample['dir'],'grouped_deletions.tsv')
        contigs = {contig.name:contig for contig in SeqIO.parse(s.references[strain],'fasta')}
        if os.path.exists(f):
            df = pd.read_csv(f,sep='\t')
            for chromosome,position,length in zip(df['chromosome'],df['position'],df['length']):
                sequence = str(contigs[chromosome][position-50:position+length+50].seq)
                #Iterating over every deletions
                if len(sequence) > 0:
                    gc_contents.append(get_gc_content(sequence))
    #Creating histogram
    fig.append_trace(go.Histogram(x=gc_contents,name='gc content',nbinsx=20),1,2)
    return fig

def main():
    """Parses deletions, gets reference gc contents and creates subplots."""
    #I dropped gc contents of the references as json to save time.
    gc_contents = get_gc_contents()

    #Creating figures for every strain
    for strain in s.strains:
        print(strain)
        #Creating subplots
        fig = make_subplots(1,2,subplot_titles=("Reference","Deletions"))
        #Adding reference gc content trace
        fig = gc_histogram(fig,gc_contents,strain)
        #Adding deletions gc content trace
        fig = deletion_histogram(fig,strain)
        fig.update_layout(
            title = strain,
            xaxis_title = 'GC content in %',
            yaxis_title = 'Counts'
        )
        fig.update_traces(showlegend=False)
        #Write to plots directory
        fig.write_image('../plots/gc_content_deletions/'+strain.replace(' ','_')+'.png')


"""
Here are some old plot functions from diversity.py that I don't use anymore but still want
to keep. In case of usage the need to be copied back do diversity.py


def ct_scatter(f, y_label):

    df = get_variants(f, platform='illumina')
    specie = 'ct'
    df = df[df['strain'] == strains[s.abbreviations[specie]]]
    out = pd.DataFrame(columns=['species', 'timepoint', 'treatment', 'snps'])
    i = 0
    for t in s.treatments[s.abbreviations[specie]]:
        for j in ['T11', 'T22', 'T33', 'T44']:
            snps = df[(df['strain'] == strains[s.abbreviations[specie]]) & (
                df['treatment'] == t) & (df['timepoint'] == j)]
            out.loc[i] = [strains[s.abbreviations[specie]],
                          j, str(t), sum(snps['hill'])]
            i += 1
    out = out.sort_values(by='treatment', ascending=True)
    fig = px.scatter(out, x='timepoint', y='snps', width=300, height=300,
                     color='treatment', color_discrete_sequence=colors['ct'])
    titles = ['T11', 'T22', 'T33', 'T44']

    for i, t in enumerate(fig['layout']['annotations']):
        t['text'] = titles[i]

    fig.for_each_yaxis(lambda y: y.update(title=''))
    fig['layout']['yaxis']['title']['text'] = y_label
    fig.update_xaxes(title='Timepoint')
    fig['layout']['legend']['title']['text'] = 'Treatment'
    fig.for_each_yaxis(lambda yaxis: yaxis.update(rangemode="tozero"))
    fig = font_size(fig)
    for d in fig['data']:
        d['marker']['size'] = 5
    if 'snps' in f:
        n = 'fixed_variants_ct_summed.svg'
    else:
        n = 'variants_ct_summed.svg'
    fig.write_image(join('..', 'plots', 'plots', n))
    return fig

def snp_distribution(f, abb, timepoint, subset=False, add_clusters=False):
    #Plots variant distributions over genome.
    #Allows to annotate "clusters". Subsetting possible with treatment
    #as string.
    hill = pd.read_csv(f, dtype={'cosm': str, 'treatment': str})
    hill = hill[hill['strain'] == s.abbreviations[abb]]
    hill = hill[hill['timepoint'] == timepoint]
    hill.index = range(len(hill))
    hill.insert(0, 'ID', hill.index)
    if subset:
        hill = hill[hill['treatment'] == subset]

    hover_data = ['pos', 'depth', 'qual', 'cosm', 'ID']

    if subset:
        n_colors = len(set(hill['cosm']))
        colors = px.colors.sample_colorscale(
            "Agsunset", [n/(n_colors - 1) for n in range(n_colors)])
        fig = px.scatter(hill, x='pos', y='freq', color='cosm', color_discrete_sequence=colors,
                         facet_col='chrom', facet_col_wrap=1,
                         facet_row_spacing=0.12, hover_data=hover_data, width=400, height=300)
    else:
        n_colors = len(set(hill['cosm']))
        colors = px.colors.sample_colorscale(
            "Agsunset", [n/(n_colors - 1) for n in range(n_colors)])
        fig = px.scatter(hill, x='pos', y='freq', color='treatment', color_discrete_sequence=colors,
                         facet_col='chrom', facet_col_wrap=1,
                         facet_row_spacing=0.12, hover_data=hover_data, width=400, height=300)
    fig.update_xaxes(matches=None, showticklabels=True)
    fig.for_each_yaxis(lambda y: y.update(title=''))
    fig['layout']['yaxis']['title']['text'] = 'Variant frequency'
    fig['layout']['xaxis']['title']['text'] = 'Position'
    for i, t in enumerate(fig['layout']['annotations']):
        t['text'] = t['text'].replace('chrom=', '')

    if add_clusters:
        # Annotations for plot
        annot = pd.read_csv(
            join('..', 'variants', 'clusters', 'clusters_'+abb+'.csv'))
        for i, chrom in enumerate(sorted(set(annot['chrom']))):
            annot_sub = annot[annot['chrom'] == chrom]
            for id in set(annot_sub['id']):
                asid = annot_sub[annot_sub['id'] == id]
                asid.index = range(len(asid))
                x0 = asid.loc[0].pos
                x1 = asid.loc[len(asid)-1].pos
                fig.add_vline(x=(x0+x1)/2, row=i, annotation_text=id,
                              opacity=0.25)
    if subset:
        n = abb+'_positions_cosm.svg'
        fig.update_layout(
            title='Variants colored by microcosm for treatment 4 in '+strains[s.abbreviations[abb]])
        fig['layout']['legend']['title']['text'] = 'Microcosm'

    else:
        fig.update_layout(title='Variant distribution along the genome')
        n = abb+'_positions.svg'
        fig['layout']['legend']['title']['text'] = 'Treatment'
    fig = font_size(fig)
    fig.for_each_yaxis(lambda yaxis: yaxis.update(
        tickmode='linear', dtick=0.2))
    fig.write_image(join('..', 'plots', 'plots',
                    n))
    return fig


def cluster_trajectories(abb):
    #This allows to follow the trajectories of SNP clusters of interest.
    df = pd.read_csv(join('..', 'variants', 'ns_variants_comp_mapping.csv'), dtype={
                     'cosm': str, 'treatment': str})
    df = df[df['strain'] == s.abbreviations[abb]]
    cluster = pd.read_csv('clusters_'+abb+'.csv')
    out = pd.DataFrame(columns=['treatment', 'chrom', 'pos', 'cosm',
                                'timepoint', 'freq', 'qual', 'linegroup', 'cluster', 'strain'])
    k = 0
    for i, c, p in zip(cluster['id'], cluster['chrom'], cluster['pos']):
        df_c = df[df['chrom'] == c]
        df_p = df_c[df_c['pos'] == p]
        for j, row in df_p.iterrows():
            out.loc[k] = [row['treatment'], row['chrom'], row['pos'],
                          row['cosm'], row['timepoint'], row['freq'], row['qual'], row['linegroup'], i, 'ct']
            k += 1
    n_colors = len(set(cluster['id']))
    if n_colors == 1:
        colors = px.colors.sample_colorscale(
            "Agsunset", [1])
    else:
        colors = px.colors.sample_colorscale(
            "Agsunset", [n/(n_colors - 1) for n in range(n_colors)])
    fig = px.line(out, x='timepoint', y='freq', line_group='linegroup',
                  facet_col='treatment', facet_col_wrap=4, color='cluster', color_discrete_sequence=colors, width=400, height=350,
                  category_orders={'timepoint': ['T11', 'T22', 'T33', 'T44'], 'treatment': sorted(set(out['treatment']))}, markers=True)
    titles = ['Treatment ' + str(i)
              for i in sorted(list(set(df['treatment'])))]
    for i, t in enumerate(fig['layout']['annotations']):
        t['text'] = titles[i]
    fig.for_each_yaxis(lambda y: y.update(title=''))
    fig['layout']['yaxis']['title']['text'] = 'Variant frequency'
    fig.update_xaxes(title=None)
    fig.for_each_xaxis(lambda x: x.update(title=''))
    fig['layout']['xaxis']['title']['text'] = 'Variant frequency'
    fig = font_size(fig)
    fig.show()


def annotate_clusters(abb):
    #Prints annotations for identified SNPs.
    def annotate_pos(gbk, c, p):
        for (start, end), (gene, product) in gbk[c].items():
            if p in range(start, end):
                return [gene, product]
        return False

    f = join('..', 'annotations', abb + '.tsv')
    df = pd.read_csv(f, sep='\t')
    # For plotting we renamed contigs to at_0 etc.
    # Rename of contigs in annotations for hashing.
    contigs = {c: abb+'_'+str(i)
               for i, c in enumerate(sorted(set(df['Sequence Id'])))}
    for i, chrom in enumerate(df['Sequence Id']):
        df.at[i, 'Sequence Id'] = contigs[chrom]
    gbk = {contig: {} for contig in sorted(set(df['Sequence Id']))}
    for i, row in df.iterrows():
        if pd.isna(row['Gene']):
            gene = 'Unknown'
        else:
            gene = row['Gene']
        if pd.isna(row['Product']):
            product = 'hypothetical protein'
        else:
            product = row['Product']
        gbk[row['Sequence Id']][(row['Start'], row['Stop'])] = (gene, product)

    out = pd.DataFrame(columns=['cluster', 'chrom', 'pos', 'gene', 'product'])
    # cluster = pd.read_csv('clusters_'+abb+'.csv')
    cluster = pd.read_csv(
        snps=join('..', 'variants', 'snps_freebayes_comp_mapping.csv'))
    for j, (i, c, p) in enumerate(zip(cluster['id'], cluster['chrom'], cluster['pos'])):
        a = annotate_pos(gbk, c, p)
        if a:
            pass
        else:
            a = ['Not annotated', 'Not annotated']
        out.loc[j] = [i, c, p, a[0], a[1]]
    out = out.drop_duplicates(subset=['cluster', 'product'])

    print(out.to_string(index=False))
"""

df = pd.read_csv(join('..','annotations','ct_variants_annotations.csv'))
df = df[df['timepoint'] == 'T44']
df = df[df['freq'] >= 0.97]
df = df[['product','treatment','cosm']]
df = df.drop_duplicates()
pac = pd.read_csv(join('..','annotations','ct_pacbio_annotations.csv'))
j = 0
k = 0
for i,row in df.iterrows():
    treatment,product,cosm = row['treatment'],row['product'],row['cosm']
    mask = (pac['treatment'] == treatment) & (pac['product'] == product) & (pac['cosm'] == cosm)
    j += 1
    if len(pac[mask]) > 0:
        k += 1
    
    