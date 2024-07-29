"""create 61x20 for subs.csv file, and creates png with the heatmap"""
import os.path

import pandas as pd
# import seaborn.apionly as sns
import matplotlib.pyplot as plt
from itertools import groupby
import matplotlib as mpl
import numpy as np
from matplotlib.patches import Rectangle


mpl.rcParams.update(mpl.rcParamsDefault)

GENCODE = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
    'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W'}
REV_GENCODE = {
    'I': ['ATA', 'ATC', 'ATT'], 'M': ['ATG'], 'T': ['ACA', 'ACC', 'ACG', 'ACT'],
    'N': ['AAC', 'AAT'], 'K': ['AAA', 'AAG'], 'S': ['AGC', 'AGT', 'TCA', 'TCC', 'TCG', 'TCT'],
    'R': ['AGA', 'AGG', 'CGA', 'CGC', 'CGG', 'CGT'], 'L': ['CTA', 'CTC', 'CTG', 'CTT', 'TTA', 'TTG'],
    'P': ['CCA', 'CCC', 'CCG', 'CCT'], 'H': ['CAC', 'CAT'], 'Q': ['CAA', 'CAG'],
    'V': ['GTA', 'GTC', 'GTG', 'GTT'], 'A': ['GCA', 'GCC', 'GCG', 'GCT'], 'D': ['GAC', 'GAT'],
    'E': ['GAA', 'GAG'], 'G': ['GGA', 'GGC', 'GGG', 'GGT'], 'F': ['TTC', 'TTT'], 'Y': ['TAC', 'TAT'],
    'C': ['TGC', 'TGT'], 'W': ['TGG']}
ALL_AA = ['A', 'Q', 'Y', 'W', 'N', 'S', 'R', 'G', 'M', 'H',
          'C', 'V', 'K', 'F', 'T', 'E', 'D', 'I/L', 'P']
ALL_CODONS = ['ATA', 'ATC', 'ATT', 'ATG', 'ACA', 'ACC', 'ACG', 'ACT', 'AAC', 'AAT', 'AAA', 'AAG',
              'AGC', 'AGT', 'AGA', 'AGG', 'CTA', 'CTC', 'CTG', 'CTT', 'CCA', 'CCC', 'CCG', 'CCT',
              'CAC', 'CAT', 'CAA', 'CAG', 'CGA', 'CGC', 'CGG', 'CGT', 'GTA', 'GTC', 'GTG', 'GTT',
              'GCA', 'GCC', 'GCG', 'GCT', 'GAC', 'GAT', 'GAA', 'GAG', 'GGA', 'GGC', 'GGG', 'GGT',
              'TCA', 'TCC', 'TCG', 'TCT', 'TTC', 'TTT', 'TTA', 'TTG', 'TAC', 'TAT', 'TGC', 'TGT',
              'TGG']  # , 'TAA', 'TAG', 'TGA']

destination_col = 'destination'
codon_col = 'codon'
protein_col = 'Leading razor protein'
position_col = 'position'
MASK_PATH = 'SixtyTwentyMask.csv'


def choose_first(stringy):
    return stringy.split(' ')[0]


def read_relevant_cols(df):
    df = df[[codon_col, destination_col, protein_col, position_col]]
    df.drop_duplicates(inplace=True)
    df = df[[codon_col, destination_col]]
    df.loc[:, destination_col] = df[destination_col].map(I_L_2IL)
    return df


def count_occurrences(mut_data, codon):
    """take a table of base and destinations codon and return for specific codon
    how many times it was mistaken for each of the amino acids"""
    codon_sliced = mut_data[destination_col][mut_data[codon_col] == codon]
    counts = codon_sliced.value_counts()
    return dict(counts)


def I_L_2IL(AA):
    if AA == 'L' or AA == 'I':
        AA = 'I/L'
    return AA


def T_2_U(codon):
    return codon.replace('T', 'U')


def U_2_T(codon):
    return codon.replace('U', 'T')


def sorting(codon):
    sam = 0
    for i in range(3):
        imp = codon[i]
        if imp == 'U' or imp == 'T': sam += 1
        if imp == 'C': sam += 2
        if imp == 'A': sam += 3
        if imp == 'G': sam += 4
        sam = sam * 10
    return sam


def add_line(ax, xpos, ypos):
    line = plt.Line2D([ypos, ypos + .2], [xpos, xpos], color='black', transform=ax.transAxes)
    line.set_clip_on(False)
    ax.add_line(line)


def label_len(my_index, level):
    labels = my_index.get_level_values(level)
    return [(k, sum(1 for i in g)) for k, g in groupby(labels)]


def label_group_bar_table(ax, df):
    xpos = -.2
    scale = 1. / df.index.size
    for level in range(df.index.nlevels):
        pos = df.index.size
        for label, rpos in label_len(df.index, level):
            add_line(ax, pos * scale, xpos)
            pos -= rpos
            lypos = (pos + .2 * rpos) * scale
            ax.text(xpos + .1, lypos, label, ha='center', transform=ax.transAxes)
        add_line(ax, pos * scale, xpos)
        xpos -= .2
    sns.set_style('darkgrid')


def create_multindex_aa_codons(df):
    index = pd.MultiIndex.from_tuples([(x, GENCODE[x]) for x in df.index],
                                      names=['codons', 'amino acid'])
    sixty_twenty = df.set_index(index)
    sixty_twenty = sixty_twenty.sort_index(level=1)
    return sixty_twenty


def plot_61x19(sixty_ninteen, plot_name, data_place, mask, nce_idx):
    fig, ax = plt.subplots(figsize=(8, 16))
    sns.set_context("poster", font_scale=0.7)
    ax = sns.heatmap(sixty_ninteen,
                     square=True,
                     cmap='Reds',
                     ax=ax,
                     cbar_kws={"shrink": .3}, linewidths=0.5,  # annot=True,
                     mask=mask)
    ax.set_facecolor('gray')
    # Below 3 lines remove default labels
    ylabels = ['' for item in ax.get_yticklabels()]
    ax.set_yticklabels(ylabels)
    ax.set_ylabel('')
    label_group_bar_table(ax, sixty_ninteen)
    fig.subplots_adjust(bottom=.001 * 2)
    for i in nce_idx:
        new_i = tuple(j+0.45 for j in i)
        ax.add_patch(Rectangle(new_i, width=0.25, height=0.25
                               , edgecolor=None
                               ,lw=0))
    plt.show()
    # plt.tight_layout()
    fig.savefig(os.path.join(data_place, plot_name + ".png"), dpi=200, bbox_inches='tight')


def ecoli_like_table(sixty_twenty):
    from matplotlib.colors import LinearSegmentedColormap
    sixty_twenty.index = list(map(T_2_U, sixty_twenty.index))
    sixty_twenty = sixty_twenty.loc[sorted(list(sixty_twenty.index), key=sorting)]
    sixty_twenty = np.log2(sixty_twenty + 1)
    plt.figure(figsize=(8, 20))

    colors = [(0, 0, 0), (1, 0, 0)]  # first color is black, last is red
    cm = LinearSegmentedColormap.from_list(
        "Custom", colors, N=20)

    sns.heatmap(sixty_twenty, cmap=cm,
                linewidths=2,
                linecolor="black",
                square=True)
    plt.savefig('61x19_ernst_design.png')


def make_sixty_twenty(mut_data):
    sixty_twenty = pd.DataFrame(columns=ALL_AA, index=ALL_CODONS)
    for codon in GENCODE:
        codon_count_dict = count_occurrences(mut_data, codon)
        for AA in codon_count_dict:
            sixty_twenty.loc[codon, AA] = codon_count_dict[AA]
    # cosmetic changes:
    sixty_twenty = sixty_twenty[sorted(list(sixty_twenty.columns),
                                       key=lambda x: 'L' if x == 'I/L' else x)]
    sixty_twenty = sixty_twenty.fillna(0)
    return sixty_twenty


def open_and_reformat_mask():
    mask = pd.read_csv(MASK_PATH, index_col=0)
    mask.index = list(map(U_2_T, mask.index))
    mask = create_multindex_aa_codons(mask)
    return mask


def nece_idx(sixty_twenty_idx, sixty_twenty_cols):
    nece = pd.read_csv("NeCEMask.csv", index_col=0)
    #sort it as the 61x20
    nece = nece.reindex(sixty_twenty_cols, axis=1)
    nece = nece.reindex(sixty_twenty_idx)
    # make a list of the index of NeCE
    nece_idx = []
    for i in range(nece.shape[0]):
        for j in range(nece.shape[1]):
            if nece.iloc[i, j]:
                nece_idx += [(j, i)]
    return nece_idx


def main(subs, output_path):
    mut_data = read_relevant_cols(subs)
    # output_path = os.path.join('..', 'results')

    sixty_twenty = make_sixty_twenty(mut_data)
    sixty_twenty = create_multindex_aa_codons(sixty_twenty)
    sixty_twenty.to_csv(os.path.join(output_path, 'sixtyTwenty.csv'))
    mask = open_and_reformat_mask()
    nece = nece_idx(sixty_twenty.index.get_level_values(0), sixty_twenty.columns)

    plot_61x19(sixty_twenty, "SixtyTwenty", output_path, mask, nece)
    log_61x19 = np.log2(sixty_twenty + 1)
    plot_61x19(log_61x19, "LogSixtyTwenty", output_path, mask, nece)


if __name__ == '__main__':
    output_dir = os.path.join('..', 'results')
    data_place = os.path.join('..', 'results', 'subs.csv')
    main(output_dir, suffix)
