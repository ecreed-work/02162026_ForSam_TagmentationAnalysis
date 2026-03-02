# TagmentationAnalysis

This Jupyter notebook analyzes tagmentation sequencing data, performing alignment, quantification, and visualization of insertion patterns. It is designed to be modular and adaptable to different experimental setups.

## Features

- Aligns reads to reference sequences using `bbmap` and `samtools`
- Extracts and quantifies integration indices
- Plots strand-specific insertion site histograms
- Outputs plots and data summaries to a specified directory

## Getting Started

### Installation

Clone this repository and create a new conda environment:

    # STEP 1: Create and activate environment
    conda create -n tagmentation python=3.10 -c conda-forge -c bioconda -y
    conda activate tagmentation
    
    # STEP 2: Install all required packages in one command
    conda install -c conda-forge -c bioconda pandas numpy matplotlib biopython logomaker openpyxl jupyterlab notebook bwa samtools
    
    # STEP 3: Install pip-only packages
    pip install pysam cutadapt
    
## Activate conda environment
Once created, you do not need to follow installation steps. Simply copy and paste:

    conda activate tagmentation

## Launch the Notebook

launch Jupyter:

    jupyter notebook

If you encounter an architecture-related error, try:

PYTHONNOUSERSITE=1 conda run -n integration_env jupyter notebook

Then, navigate to the location of the notebook and open CAST_TagmentationNotebook.ipynb.

## Usage
First use map_insertions_template.xlsx
    replace path to fastq files (these come from MiSeq)
    Change reference genome file path
    indicate path to directory you would like tsv files to be output to

Start with map_insertions.py
    at top of script, there is a commented out section
    this section contains commands you can edit, copy, and paste into the terminal
    replace path to xlsx file with your own, tells the code how to find fastq files from miseq and where to place output tsv files

Once insertions have been mapped, tsv_analysis.ipynb
    Select bulk or small scale input
    run cells to generate insertion profiles based on bulk of small scale use, labeled at top
    after genome insertion profiles are generated, run remaining cells sequentially

To analyze 

1. Modify the input parameters (e.g., input FASTA/FASTQ files, output directory) where comments specify in the notebook.

2. Run each cell sequentially.

3. Output files including histograms and processed data will be saved in the specified output directory.

## Dependencies

This notebook requires the following packages (automatically installed with the conda environment above):

biopython
matplotlib
numpy
bbmap
samtools
notebook
jupyterlab
