# -*- coding: utf-8 -*-
"""
Created on 2024

@author: ernestmordret, ommerasraf, tehilaleiman

"""

import numpy as np
import pandas as pd
from pathlib import Path
import subsUtils
from Bio import SeqIO
from Bio import Seq


def parse_list_in_df(df_cell):
    return eval('[' + df_cell.replace(';', ',') + ']')


def get_column_names(dp_path):
    """because the name of the columns can be changed according to the samples names"""
    dp_cols = pd.read_csv(dp_path, sep='\t', nrows=0)
    bp_cols = list(filter(lambda x: 'Base intensity' in x, dp_cols))
    dp_cols = ['Localization probability', 'Modification intensities', 'Base peptide sequence', 'Modification',
               'Cluster mass', 'Mass difference', 'Time difference', 'Posterior error probability', 'Score',
               'Max Probability'] + bp_cols
    return dp_cols, bp_cols


def split_clusters(dp):
    # clusters in dependentPeptides table are defined by base peptide sequence and modification
    # a problem occurs when the modification is 'Unmodified' and then peptides with different mass differences are
    # clustered together. I don't trust it fully, so I'm taking apart the clusters
    dp['Mass difference'] = dp['Mass difference'].map(parse_list_in_df)
    dp['Time difference'] = dp['Time difference'].map(parse_list_in_df)
    dp['Modification intensities'] = dp['Modification intensities'].map(parse_list_in_df)
    dp = dp.explode(['Mass difference', 'Modification intensities', 'Time difference'])
    dp['Mass difference'] = dp['Mass difference'].astype(float)
    dp['Time difference'] = dp['Time difference'].astype(float)
    dp['dp_intensity'] = dp['Modification intensities'].astype(float)
    return dp


def filter_localization_threshold(dp):
    # filter the AAs that were possibly modified (with p > threshold)
    dp = dp[dp['Max Probability'] >= subsUtils.AA_LOC_THRESH]
    dp.loc[:, ['error_position_in_pep', 'origin']] = dp['Localization probability'].map(
        subsUtils.refine_localization_probabilities).to_list()
    dp = dp[dp['error_position_in_pep'].notna()]
    return dp


def map_to_prot_and_filter_contaminants(path_to_peptides, dp):
    peptides_cols = ["Reverse", "Potential contaminant", "Sequence", "Start position", "Leading razor protein",
                     "Unique (Proteins)"]

    # filter out contaminants and decoys, and add information
    peptides = pd.read_csv(path_to_peptides, sep="\t", usecols=peptides_cols)
    dp = dp.merge(peptides, how="left", left_on="Base peptide sequence", right_on="Sequence")
    dp = dp[dp["Potential contaminant"] != "+"]
    dp = dp[dp["Reverse"] != "+"]  # not suppose to exist, but better safe than sorry
    dp = dp[~(dp["Leading razor protein"].str.contains(
        'REV'))]  # for some reason, the filter for reverse itself is not good enough
    dp['Start position'] = dp['Start position'].map(int) - 1
    return dp


def remove_duplicates_from_fasta(fasta_file):
    """This function reads a FASTA file, removes duplicate gene sequences based on both the gene name (identifier) and the sequence,
    and writes the unique sequences to a new output file. Sequences are considered duplicates if both the identifier and the sequence
    are identical.
    """
    unique_sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        identifier = str(record.id)
        sequence = str(record.seq).upper()

        # Use a combination of the identifier and sequence as the key to avoid duplicates
        if (identifier, sequence) not in unique_sequences:
            unique_sequences[(identifier, sequence)] = record
    return unique_sequences.values()


#
# def remove_duplicates_from_fasta(fasta_file):
#     """This function reads a FASTA file, removes duplicate gene sequences based on the gene identifier,
#     and writes the unique sequences to a new output file. If there are multiple sequences with the same identifier
#     but different sequences, only the first occurrence is retained."""
#     unique_identifiers = set()
#     unique_sequences = {}
#
#     for record in SeqIO.parse(fasta_file, "fasta"):
#         identifier = str(record.id)
#
#         # Check if the identifier is already processed
#         if identifier not in unique_identifiers:
#             unique_identifiers.add(identifier)
#             # Use the identifier as the key to avoid duplicates based on it
#             unique_sequences[identifier] = record
#     return unique_sequences.values()


def map_protein_to_sequence(dna_fasta_path, aa_fasta_path, dp):
    dna_mapped, aa_mapped, dna_str_fasta, aa_str_fasta = pd.DataFrame(), pd.DataFrame(), '', ''
    len_dp = len(dp)
    dp["mapped"] = False
    dp["dna_seq"] = None
    if dna_fasta_path:
        fasta = remove_duplicates_from_fasta(dna_fasta_path)
        dna_records = SeqIO.to_dict(fasta)
        dp["mapped"] = dp["Leading razor protein"].map(lambda x: x in dna_records)  # only proteins that exist in fasta
        dna_mapped = dp[dp["mapped"]].copy()
        dna_mapped["dna_seq"] = dna_mapped["Leading razor protein"].map(lambda geneid: dna_records[geneid].seq)
        dna_mapped["protein_seq"] = dna_mapped["dna_seq"].map(lambda seq: seq.translate() if seq else None)
        dna_str_fasta = ''.join(map(lambda geneid: str(dna_records[geneid].seq.translate), dna_records))
    dp_unmapped = dp[~dp["mapped"]]
    if aa_fasta_path and not dp_unmapped.empty:
        fasta = remove_duplicates_from_fasta(aa_fasta_path)
        aa_records = SeqIO.to_dict(fasta)
        dp_unmapped["mapped"] = dp_unmapped["Leading razor protein"].map(
            lambda x: x in aa_records)  # only proteins that exist in fasta
        aa_mapped = dp_unmapped[dp_unmapped["mapped"]].copy()
        aa_mapped["protein_seq"] = aa_mapped["Leading razor protein"].map(lambda geneid: aa_records[geneid].seq)
        dna_str_fasta = ''.join(map(lambda geneid: str(aa_records[geneid].seq) + '*', aa_records))

    dp_mapped = pd.concat([dna_mapped, aa_mapped])
    print(f'{len(dp_mapped)} were mapped to fasta file out of {len_dp} rows '
          f'({(len(dp_mapped) / len_dp) * 100: .1f}%)')
    str_fasta = dna_str_fasta + aa_str_fasta
    return dp_mapped, str_fasta


def get_error_position(dp):
    # Does the peptide come from the N/C-term of the protein
    dp["position"] = (dp["Start position"] + dp["error_position_in_pep"]).map(int)
    dp['prot_terminal'] = dp.apply(
        lambda row: subsUtils.peptide_termini(len(row['protein_seq']), row['error_position_in_pep']), axis=1)
    return dp


def get_known_ptms(dp, tol):
    # Handles mass differences that can be explained as PTMs
    danger_mods = pd.read_pickle('danger_mods_new_unimods.pkl')
    dp['ptm'] = False
    for ind, mod in danger_mods.iterrows():
        position = mod['position']
        site = mod['site']
        delta_m = mod['delta_m']

        # tolerance of mass here is twice than filtering for aa mass change
        mass_filter = (delta_m - (2 * tol) < dp['Mass difference']) & \
                      (dp['Mass difference'] < delta_m + (2 * tol))

        term_filters_dict = {'Protein N-term': dp['prot_terminal'] == 'n-term',
                             'Protein C-term': dp['prot_terminal'] == 'c-term',
                             'Any N-term': dp['pep_terminal'] == 'n-term',
                             'Any C-term': dp['pep_terminal'] == 'c-term',
                             'Anywhere': True}
        term_filter = term_filters_dict[position]

        site_filter = True
        if site in subsUtils.MW_DICT.keys():
            site_filter = dp["origin"] == site

        dp.loc[site_filter & term_filter & mass_filter, 'ptm'] = True
    dp = dp[~dp['ptm']]
    return dp


def get_dp_table_and_clean(dp_path, pep_path, dna_fasta_path, aa_fasta_path, tol):
    dp_cols, bp_cols = get_column_names(dp_path)
    dp = pd.read_csv(dp_path, sep='\t', usecols=dp_cols)
    dp = split_clusters(dp)
    dp = filter_localization_threshold(dp)
    dp = map_to_prot_and_filter_contaminants(pep_path, dp)

    # annotate whether the error was in the peptide termini
    dp['peptide_length'] = dp['Base peptide sequence'].map(len)
    dp['pep_terminal'] = dp.apply(
        lambda row: subsUtils.peptide_termini(row['peptide_length'], row['error_position_in_pep']), axis=1)

    dp, str_fasta = map_protein_to_sequence(dna_fasta_path, aa_fasta_path, dp)
    dp = get_error_position(dp)

    # Handles mass differences that can be explained as PTMs
    dp = get_known_ptms(dp, tol)
    return dp, bp_cols, str_fasta


def get_masses_that_match_substitutions(dp, tol):
    dp['substitution'] = False
    subs_dict = subsUtils.get_subs_dict()
    for i in sorted(subs_dict.keys()):
        delta_m = subs_dict[i]
        original_aa = i[0]
        # TEHILA: for each mutation slice the substitutions that have happened to specific amino acid,
        # and have the expected mass difference for the mutation, and the probability that it is a dependant peptide
        # is 0.95.
        dp.loc[
            (dp['Mass difference'] > delta_m - tol) &
            (dp['Mass difference'] < delta_m + tol) &
            (dp['origin'] == original_aa),
            'substitution'] = i
    subs = dp[dp['substitution'] != False].copy()
    return subs


def map_subs(subs):
    subs['codon'] = subs.apply(lambda x: x['dna_seq'][
                                         x['position'] * 3:(x['position'] * 3) + 3
                                         ] if isinstance(x['dna_seq'], Seq.Seq) else None, axis=1)
    subs['destination'] = subs['substitution'].map(lambda x: x[-1] if x else False)

    mask = subsUtils.get_mispairing_mask()
    subs['mispairing'] = subs.apply(subsUtils.is_mispairing, axis=1, args=(mask,))

    subs['dp_sequence'] = subs.apply(lambda x: subsUtils.create_modified_seq(x['Base peptide sequence'],
                                                                             x['error_position_in_pep'],
                                                                             x['destination']),
                                     axis=1)
    return subs


def handling_leucine_isoleucine(subs, str_fasta):
    subsl = subs[subs['destination'] == 'L']
    subsnol = subs[subs['destination'] != 'L']
    subsl['dp_sequence'] = subsl.apply(lambda x: subsUtils.create_modified_seq(x['Base peptide sequence'],
                                                                               x['error_position_in_pep'],
                                                                               x['destination'],
                                                                               'I'),
                                       axis=1)
    subsl = subsl[subsl['dp_sequence'].map(lambda seq: seq not in str_fasta)]
    subs = pd.concat([subsnol, subsl], ignore_index=True)
    return subs


def calculate_dp_intensity(subs):
    grp_by = ['Sequence', 'dp_sequence', 'origin', 'destination', 'substitution', 'mispairing', 'Leading razor protein',
              'position', 'error_position_in_pep', 'Start position', 'peptide_length', 'Posterior error probability',
              'Score', 'codon', 'Localization probability', 'Max Probability', 'Modification', 'dna_seq',
              'protein_seq', 'bp_intensity', 'Unique (Proteins)']
    subs = subs.replace(np.nan, 'None')  # otherwise it drops with groupby

    # subs["Modification"] = subs.Modification.replace(np.nan, 'None')  # otherwise it drops with groupby
    subs = subs.groupby(grp_by).agg({'Mass difference': 'mean', 'Time difference': 'mean',
                                     'dp_intensity': 'sum'}).reset_index()
    subs = subs.rename(columns={'Sequence': 'bp_sequence', 'Start position': 'peptide_start_position'})

    # After the grouping step, two types of duplicates can still persist:
    # 1. Cases where the modification has different PEP and localization probability.
    # 2. Instances where the same mass of DP is identified by MS but is attributed to two different BPs. In such cases,
    # the same modified sequence is detected by our pipeline.
    # To address these duplicates, I opt to retain only one instance of each DP. I choose to keep the one with the
    # higher BP intensity and higher PEP, assuming it as a reasonable strategy.

    subs = subs.sort_values(by=['bp_intensity', 'Posterior error probability']).drop_duplicates('dp_sequence')
    return subs


def get_intensities_values(subs, bp_cols):
    # the mean isn't really meaningful because they all have the same value, so instead of choosing one of them the
    # code is calculating the mean - which is the same value
    subs['bp_intensity'] = subs[bp_cols].mean(axis=1)
    subs = subs.drop(columns=bp_cols)
    subs = calculate_dp_intensity(subs)
    # add ratio col
    subs['dp/bp_intensities_ratio'] = subs['dp_intensity'] / subs['bp_intensity']
    subs['log_intensities_ratio'] = np.log10(subs['dp/bp_intensities_ratio'])
    return subs


def clean_and_process_subs(subs, str_fasta, bp_cols):
    subs = map_subs(subs)
    # remove peptides that exist in the protein fasta
    subs = subs[subs['dp_sequence'].map(lambda seq: seq not in str_fasta)]
    subs = handling_leucine_isoleucine(subs, str_fasta)
    subs = subs[subs['Posterior error probability'] <= 0.05]
    subs = get_intensities_values(subs, bp_cols)
    return subs


def save_subs(subs, output_dir):
    if not output_dir:
        output_dir = Path('output')
    output_dir = Path(output_dir)
    if output_dir.exists():
        print('DETECT:    output dir exists, overwriting')
    output_dir.mkdir(exist_ok=True)
    # subs.to_pickle(output_dir / 'subs')
    subs.to_csv(output_dir / 'subs.csv')


# %%
def main(dependent_peptides_path, peptides_path, dna_fasta_path=None, aa_fasta_path=None, output_dir=None):
    if (not dna_fasta_path) and (not aa_fasta_path):
        print("No FASTA file was provided. To run this code, a FASTA file with headers identical to the ones that"
              "has been used by MaxQuant should be provided.")
        return
    current_working_directory = Path.cwd()
    print("The location of the current working directory is:\n", current_working_directory)

    # Initialize variables
    path_to_dependentpeptides = Path(dependent_peptides_path)
    path_to_peptides = Path(peptides_path)

    # Define constants
    tol = 0.005  # tolerance for variation of mass

    # %% code
    dp, bp_cols, str_fasta = get_dp_table_and_clean(path_to_dependentpeptides, path_to_peptides,
                                                    dna_fasta_path, aa_fasta_path, tol)
    subs = get_masses_that_match_substitutions(dp, tol)
    if subs.empty:
        print("No errors were found")
        return
    subs = clean_and_process_subs(subs, str_fasta, bp_cols)
    save_subs(subs, output_dir)
    return subs


if __name__ == '__main__':
    DATA_PATH = Path(r'../data')
    DP_PATH = DATA_PATH / 'dependentPeptides.txt'
    PEP_PATH = DATA_PATH / 'peptides.txt'
    DNA_FASTA_PATH = DATA_PATH / 'Homo_sapiens.GRCh38.cds.all.fa'
    AA_FASTA_PATH = DATA_PATH / 'new.modheader.all_nonsyn.fasta'
    OUTPUT_PATH = Path(r'../results')

    subs = main(dependent_peptides_path=DP_PATH,
                peptides_path=PEP_PATH,
                dna_fasta_path=DNA_FASTA_PATH,
                aa_fasta_path=AA_FASTA_PATH,
                output_dir=OUTPUT_PATH)
