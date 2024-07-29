# Running MaxQuant for Input Preparation

This guide provides instructions on how to run MaxQuant to generate the necessary input files (`dependentPeptides.txt`, `peptides.txt`) for use with the LC-MS Phenotypic Variation Detector pipeline.

## Installation

We used MaxQuant 2.1.3.0

1. Download MaxQuant from [MaxQuant official website](https://cox-labs.github.io/coxdocs/Download_Installation.html).
2. Follow the installation instructions provided on the site to install MaxQuant on your system.

## Running MaxQuant

1. Open MaxQuant.
2. Configure your experiment settings, including the paths to your raw data files.
3. In the 'Global parameters' tab, choose 'Identification' and check 'Dependent peptides' box.
Note that the pipeline work best with trypsin digested proteomics data.

## Running the Analysis

If you run it on windows, click 'Start' to begin the analysis.
To run it on linux system, read the instructions on [MaxQuant official website](https://cox-labs.github.io/coxdocs/Download_Installation.html).
Once the run is complete, check the output directory 'combined\txt' for the `dependentPeptides.txt` and `peptides.txt` files. These files will be used as input for the LC-MS Phenotypic Variation Detector.
