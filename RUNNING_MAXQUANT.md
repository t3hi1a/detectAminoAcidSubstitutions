# Running MaxQuant for Input Preparation

This guide provides instructions on how to run MaxQuant to generate the necessary input files (`dependentPeptides.txt`, `peptides.txt`) for use with the LC-MS Phenotypic Variation Detector pipeline.

## Installation

We used MaxQuant 2.1.3.0

1. Download MaxQuant from [MaxQuant official website](https://cox-labs.github.io/coxdocs/Download_Installation.html).
2. Follow the installation instructions provided on the site to install MaxQuant on your system.

## Running MaxQuant

1. Open MaxQuant GUI.
2.	Works on *.raw files (note: can occupy a lot of disk space).
3.	First set the experiment: set experiment (e.g. ctrl vs treatment). If there are fractions, set fraction number for each raw file.
4.	Choose referrence FASTA at `Sequences` at `Global parameters`. There is a possibility to use the protein sequences, or the nucleotide (CDS sequences). The second is preferred, because in later stages the pipeline is able to identify the codon with the change.
5.	In the `Global parameters` tab, choose `Identification` and check `Dependent peptides` box.
6.	Other parameters can be set in consultation with an MS expert.
7.	After setting the parameters it is also possible to save the parameters to file. (file-> save parameters). This generates *.xml file which can be migrated to different server. If it is moved, the paths to the files inside the parameters file should be change. This can be done by command line MaxQuant --changeFolder

Note that the pipeline work best with trypsin digested proteomics data.

## Running the Analysis

If you run it on windows, click 'Start' to begin the analysis.
To run it on linux system, read the instructions on [MaxQuant official website](https://cox-labs.github.io/coxdocs/Download_Installation.html).
Once the run is complete, check the output directory `combined\txt` for the `dependentPeptides.txt` and `peptides.txt` files. These files will be used as input for the LC-MS Phenotypic Variation Detector.
