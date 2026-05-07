This repository contains bioinformatics analysis scripts and supporting files for the manuscript:

**“A high-penetrance intergenic variant at 9p21 confers melanoma susceptibility.”**

This study was conducted as part of the MelaNostrum Consortium:  
https://dceg.cancer.gov/research/cancer-types/melanoma/melanostrum

The repository is organized into two main analysis modules:

## 1. AlphaGenome prediction

This folder contains code and input sequence files used to predict the regulatory effects of the 9p21 deletion using the deep learning model AlphaGenome.

Included files:

- `AlphaGenome_prediction_9p21deletion.ipynb`  
  Jupyter Notebook used to run AlphaGenome predictions for the reference and deletion alleles.

- `AlphaGenome_prediction_9p21deletion.html`  
  Static rendered HTML version of the notebook for convenient viewing without running Jupyter.

- FASTA files containing GRCh38/hg38 genomic sequences used as input for AlphaGenome prediction.

## 2. Structural variant calling

This folder contains the script used to call structural variants from the whole-genome sequencing BAM file using ClinSV.

Included files:

- R script for plotting read coverage for the sample `7646-01005` 
- Bash script for running ClinSV structural variant calling.
- Input files required for the analysis.
- Information for accessing the WGS BAM file through dbGaP.

The WGS data generated for sample `7646-01005` are available through dbGaP under accession:`phs004649.v1`

Individual-level sequencing data are not included in this repository because of participant privacy and controlled-access data restrictions.

## Software dependencies

The analyses were performed using the following software:

- Python 3.10
- R 4.0
- Jupyter Notebook or JupyterLab
- AlphaGenome v0.6.1. Additional Python package dependencies for the AlphaGenome analysis are listed within the Jupyter Notebook.
- ClinSV v1.1
- Singularity v4.3.7
- Additional R package dependencies for generating plots are listed within the R script.
