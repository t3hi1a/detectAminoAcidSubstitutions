import sys
import os
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

from scripts import detect_main, confusion_matrix
import configparser

CONFIG_FILE = 'config.ini'


def load_config(config_file):
    # Load configuration file
    config = configparser.ConfigParser()
    config.read(config_file)

    # Retrieve paths from the configuration file
    output_dir = Path(config['Paths']['output_dir'])

    # Retrieve file names from the configuration file
    dp_tb = Path(config['Paths']['dependentpeptides_file'])
    pep_tb = Path(config['Paths']['peptides_file'])
    dna_file = Path(config['Paths']['dna_file'])
    proteins_file = Path(config['Paths']['proteins_file'])
    return dp_tb, pep_tb, dna_file, proteins_file, output_dir


if __name__ == '__main__':
    dp_tb, pep_tb, dna_file, proteins_file, output_dir = load_config(CONFIG_FILE)

    subs = detect_main(dp_tb, pep_tb, dna_file, proteins_file, output_dir)

    # if dna fasta was provided and the substitutions were mapped to
    subs_codons = subs[subs.codon != 'None']
    if not subs_codons.empty:
        confusion_matrix(subs_codons, output_dir)

    plt.figure()
    sns.histplot(subs['log_intensities_ratio'])
    plt.savefig(os.path.join(output_dir, "intensity_ratio_distribution.png"))