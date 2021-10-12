from samples import Samples
import subprocess

def plot_gc_bias():
    """Plotting GC bias, plotting is done by package gc_bias."""
    s = Samples()

    for strain in s.strains:
        for sample in s.strains[strain]:
            #Creating labels json for every sample
            labels = dict()
            labels['title'] = strain+' in '+sample['name']
            labels['xlabel'] = 'GC content per window'
            labels['ylabel'] = 'coverage per window'
            labels['theme'] = 'plotly_dark'

            j = json.dumps(labels,indent=4) 
            with open(os.path.join(sample['dir'],'labels.json'),'w') as handle:
                handle.write(j)

            #Calling bash script to call gc_bias package
            #See https://github.com/nahanoo/gc_bias
            cmd = ['sbatch','plot_gc_bias.sh',sample['dir']]
            print(sample['dir'])
            subprocess.call(' '.join(cmd),shell=True)
    
if __name__ == "__main__":
    plot_gc_bias()