# Project Title

# LC-MS Phenotypic Variation Detector

This project is designed to process LC-MS results, specifically those analyzed by MaxQuant, to detect amino acid substitutions within proteins and identify phenotypic variations potentially caused by mutations, transcription errors, translation errors, or post-translational modifications (PTMs). It is particularly useful for researchers focused on understanding and analyzing phenotypic variations and their underlying causes.
This project builds upon the methodologies described in the [2019 publication on PubMed](https://pubmed.ncbi.nlm.nih.gov/31353208/).

## Key Features

- **Phenotypic Variation Detection**: Automatically detects and identifies amino acid substitutions in protein sequences.
- **Comprehensive Results**: Outputs a detailed table listing all detected phenotypic variations along with intensity ratios, providing a quantitative measure of each variation.
- **Enhanced Performance**: Offers improved speed and accuracy compared to previous implementations, making it a valuable tool for extensive datasets.

## Technological Highlights

The tool is built using Python, leveraging libraries such as Pandas for data manipulation, Matplotlib for generating visual representations of the data, and BioPython for handling biological computations efficiently. This combination of technologies ensures a robust solution capable of handling complex biological data.

## Audience

This tool is intended for scientific researchers and bioinformaticians interested in exploring and understanding phenotypic variations in proteins due to various biological processes.




## Files Description

- **`danger_mods_new_unimods.pkl`**: Contains a table of protein modifications sourced from [UniMod](https://www.unimod.org/). This file helps eliminate mass changes due to post-translational modifications (PTMs) or artifacts, which are not from substitutions.

- **`detect.py`**: Processes 'dependentPeptides.txt' and 'peptides.txt' from MaxQuant output, along with a FASTA file (ideally DNA FASTA to facilitate codon identification). It outputs a table 'subs.csv' listing all detected substitutions, ensuring gene names match those used in MaxQuant.

- **`main.py`**: The central script that orchestrates the project workflow. It integrates functionalities from all other scripts, executes `detect.py`, and generates both a histogram of intensity ratios and a heatmap of substitution counts.

- **`NeCEMask.csv`** and **`SixtyTwentyMask.csv`**:
  - `NeCEMask.csv`: Identifies near-cognate substitutions where only one nucleotide change could explain the amino acid substitution.
  - `SixtyTwentyMask.csv`: Marks invalid cells either as cognate codons or indistinguishable substitutions due to mass similarities with other substitutions or modifications.

  Both CSV files are crucial for plotting the heatmap of substitutions in `subsToSixtyTwentyWithMask.py`.

- **`plot_intensity_ratio.py`**: Generates histograms of intensity ratios, illustrating the quantitative difference between peptides with and without substitutions, thereby providing insights into the impact of these variations.

- **`subsToSixtyTwentyWithMask.py`**: Processes 'subs.csv' to produce a heatmap detailing the frequency of each type of substitution (e.g., codon to amino acid) in the sample, providing a visual representation of substitution patterns.

- **`subsUtils.py`**: Contains essential utility functions that support operations in `detect.py`, including various data manipulation and processing tasks.


## Installation

Describe how to install the project, e.g.,

```bash
git clone https://yourprojectlink.git
cd yourprojectdirectory
pip install -r requirements.txt
