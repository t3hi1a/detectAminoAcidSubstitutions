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

- `danger_mods_new_unimods.pkl`: A pickled data file containing [describe the data and its usage in the project].
- `detect.py`: Script for detecting [describe what it detects and how it's used].
- `main.py`: The main entry point for the project. This script [describe its role].
- `NeCEMask.csv`: CSV file containing [describe the contents and how they are used].
- `plot_intensity_ratio.py`: This script generates plots for [describe what it plots].
- `run`: Bash script to run the project. This script simplifies the execution process.
- `SixtyTwentyMask.csv`: Contains [describe data].
- `subsToSixtyTwentyWithMask.py`: Script to [describe its function].
- `subsUtils.py`: Utility functions used across the project. It includes functions to [describe types of utilities].

## Installation

Describe how to install the project, e.g.,

```bash
git clone https://yourprojectlink.git
cd yourprojectdirectory
pip install -r requirements.txt
