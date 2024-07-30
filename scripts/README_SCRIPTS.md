## Files Description

- **`detect.py`**: Processes 'dependentPeptides.txt' and 'peptides.txt' from MaxQuant output, along with a FASTA file (ideally DNA FASTA to facilitate codon identification). It outputs a table 'subs.csv' listing all detected substitutions, ensuring gene names match those used in MaxQuant.

- **`plot_intensity_ratio.py`**: Generates histograms of intensity ratios, illustrating the quantitative difference between peptides with and without substitutions, thereby providing insights into the impact of these variations.

- **`subsToSixtyTwentyWithMask.py`**: Processes 'subs.csv' to produce a heatmap detailing the frequency of each type of substitution (e.g., codon to amino acid) in the sample, providing a visual representation of substitution patterns.

- **`subsUtils.py`**: Contains essential utility functions that support operations in `detect.py`, including various data manipulation and processing tasks.
