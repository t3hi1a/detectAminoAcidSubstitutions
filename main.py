import sys
import os
import matplotlib.pyplot as plt
import seaborn as sns
import detect
import subsToSixtyTwentyWithMask
from pathlib import Path

if __name__ == '__main__':
    #
    slava_dir = Path(r'C:\Users\tehilal\OneDrive - weizmann.ac.il\Weizmann\Pilpelab\Collateral Vaccine Design\SixtyTwenty collection\Slava')
    dp_tb = slava_dir / 'dependentPeptides.txt'
    pep_tb = slava_dir / 'peptides.txt'
    dna_file = ''
    proteins_file = slava_dir / 'proteins.fa'
    output_dir = ''


    subs = detect.main(dp_tb, pep_tb, dna_file, proteins_file, output_dir)
    subs_codons = subs[subs.codon != 'None']
    if not subs_codons.empty:
        subsToSixtyTwentyWithMask.main(subs_codons, output_dir)
    plt.figure()
    sns.histplot(subs['log_intensities_ratio'])
    plt.savefig(os.path.join(output_dir, "intensity_ratio_distribution.png"))
