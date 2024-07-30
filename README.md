# LC-MS Phenotypic Variation Detector

This project is designed to process LC-MS results, specifically those analyzed by MaxQuant, to detect amino acid substitutions within proteins, potentially caused by mutations, transcription errors, translation errors, or post-translational modifications (PTMs).
This project builds upon the methodologies described in the [2019 publication on PubMed](https://pubmed.ncbi.nlm.nih.gov/31353208/).

## Key Features

- **Phenotypic Variation Detection**: Automatically detects and identifies amino acid substitutions in protein sequences.
- **Comprehensive Results**: Outputs a detailed table listing all detected phenotypic variations along with intensity ratios, providing a quantitative measure of each variation.
- **Enhanced Performance**: Offers improved speed and accuracy compared to previous implementations, making it a valuable tool for extensive datasets.

## Audience

This tool is intended for scientific researchers and bioinformaticians interested in exploring and understanding substitution of amino acids due to various biological processes.


## Installation
first, you need to:
```bash
module load git-lfs
```
or
```bash
git lfs install
```
and then
```bash
git lfs clone https://github.com/t3hi1a/detectAminoAcidSubstitutions.git
cd detectAminoAcidSubstitutions
pip install -r requirements.txt
