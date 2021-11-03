# Sequencing data analysis for Piccardis evolution experiment

The name of this repository is depreciated and will change soon. 

This repository hosts code for analyzing the sequencing data derived from a community evolution experiment.
Four strains were cultured together over 44 weeks. As a control, if possible, the strains were also grown in monoculture.
From this experiment we have metagenomic Illumina sequencing data available for the timepoints 11, 22, 33 and 44 (in weeks).
At week 44, samples were plated and picked colonies were amplified. The DNA was isolated and sequenced with PacBio.

If you want to have a glims at the results you can do that by having a look at this [report](https://github.com/nahanoo/black_queen_hypothesis/blob/main/reports/report.pdf).  
For sample processing, a Snakemake workflow was implemented for the PacBio and the Illumina data.
The PacBio workflow can be found [here](https://github.com/nahanoo/black_queen_hypothesis/tree/main/scripts/workflows/pacbio) and the
Illumina workflow [here](https://github.com/nahanoo/black_queen_hypothesis/tree/main/scripts/workflows/illumina)

The data resulting from those workflows was analyzed in [analyze_pacbio.py](https://github.com/nahanoo/black_queen_hypothesis/blob/main/scripts/analyze_pacbio.py) and in [analyze_illumina.py](https://github.com/nahanoo/black_queen_hypothesis/blob/main/scripts/analyze_illumina.py).


For this analysis I made use of package that I wrote called [deletion_detection](https://github.com/nahanoo/deletion_detection).
It allows to identify regions with no coverage and reads with deletions. Those regions are then annotated by using a parsed
genbank file.  
Alignment files derived from Illumina data were analyzed using [gc_bias](https://github.com/nahanoo/gc_bias).

If someone is ever going to read this code you will see that I import a samples class. This is a sample parser package that I wrote.
It's very handy but contains a lot of server paths which is why it's located in a private repository.