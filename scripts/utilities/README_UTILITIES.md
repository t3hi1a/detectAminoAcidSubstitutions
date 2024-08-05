**file description**

- **`danger_mods_new_unimods.pkl`**: Contains a table of protein modifications sourced from [UniMod](https://www.unimod.org/). This file helps eliminate mass changes due to post-translational modifications (PTMs) or artifacts, which are not from substitutions.

- **`NeCEMask.csv`** and **`SixtyTwentyMask.csv`**:
  - `NeCEMask.csv`: Identifies near-cognate substitutions where only one nucleotide change could explain the amino acid substitution.
  - `SixtyTwentyMask.csv`: Marks invalid cells either as cognate codons or indistinguishable substitutions due to mass similarities with other substitutions or modifications.

  Both CSV files are crucial for plotting the heatmap of substitutions in `subsToSixtyTwentyWithMask.py`.
