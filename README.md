# Sequencing data analysis for the paper "The evolution of reduced facilitation in a four-species bacterial community"

This repository hosts the code for the sequencing data analysis.
The experiment described in the paper was followed using DNA short-read, DNA long-read and RNA short-read sequencing data.

## Environment

The conda environment containing most of the software used for data analysis is exported in `environment.yml`

## Data processing

The genomic short-read and long-read data was processed using a Snakemake workflow.  
They can be found in the following directories:
- Illumina: `scripts/workflows/illumina`
- PacBio: `scripts/workflows/pacbio`

The RNA sequencing data was analyzed usign [RASflow](https://github.com/zhxiaokang/RASflow) workflow.

## Data analysis

Variants and other information relevant for the analysis was parsed and stored ad dataframes in the `variants` directory.  
Below is an outline of the most important dataframes:
- `variants/variants_comp_mapping.csv` contains all filtered variants for the genomic Illumina data
- `variants/{at,ct}_variants_annotations.csv` same as above, but additionally with genomic annotations for Ct and At
- `variants/snps_freebayes_comp_mapping.csv` contains all fixated variants for the genomic Illumina data
- `variants/snps_pacbio.csv` contains all SNPS in the genomic PacBio data
- `variants/assembly_length.csv` assembly lenghts for PacBio assemblies

The code for the figures can be found in `scripts/diversity.py` for variant analsysis and `scripts/deletions.py` for PacBio data.
