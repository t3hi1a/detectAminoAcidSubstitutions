#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  7 18:48:17 2019

@author: omer
"""

import pandas as pd
import numpy as np
import re
from Bio import SeqIO
from itertools import groupby
from operator import itemgetter


# %% constants
MW_DICT = {
    "G": 57.02147, "A": 71.03712, "S": 87.03203, "P": 97.05277, "V": 99.06842,
    "T": 101.04768, "I": 113.08407, "L": 113.08407, "N": 114.04293, "D": 115.02695,
    "Q": 128.05858, "K": 128.09497, "E": 129.0426, "M": 131.04049, "H": 137.05891,
    "F": 147.06842, "R": 156.10112, "C": 160.030654, "Y": 163.0633, "W": 186.07932
}

PEP_FASTA = "Homo_sapiens.GRCh38.pep.all.fa"
CODONS_DICT = {
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
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
    'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W'}
CODONS = list(CODONS_DICT.keys())
AMINO_ACIDS = [CODONS_DICT[codon] for codon in CODONS]
AA_LOC_THRESH = 0.95


# %%

def get_inverted_codon_table():
    ct = CODONS_DICT
    inv_codon_table = {}
    for k, v in ct.copy().items():
        inv_codon_table[v] = inv_codon_table.get(v, [])
        inv_codon_table[v].append(k)
    inv_codon_table['L'] += inv_codon_table['I']
    return inv_codon_table


def get_subs_dict():
    subs_dict = {}
    for i in MW_DICT:
        for j in MW_DICT:
            if i != j and j != 'I':
                subs_dict[i + ' to ' + j] = MW_DICT[j] - MW_DICT[i]
    del subs_dict['I to L']
    subs_dict = {key.replace('to L', 'to I/L'): subs_dict[key] for key in subs_dict}
    return subs_dict

def codonify(seq):
    """
    input: a nucleotide sequence (not necessarily a string)
    output: a list of codons
    """
    seq = str(seq)
    l = len(seq)
    return [seq[i:i + 3] for i in range(0, l, 3)]


def find_proteins(base_seq, names_list, boundaries_aa, W_aa, sa):
    """
    input: a peptide sequence (string)
    output: the names of proteins containing that sequence 
    """
    tbr = " ".join([names_list[i] for i in np.searchsorted(boundaries_aa - 1, SA_search(base_seq, W_aa, sa))])
    if tbr.strip(" ") == '':
        return ''
    else:
        return tbr


def fetch_codon(base_seq, modified_pos, record_dict, names_list, boundaries_aa, W_aa, sa):
    """
    input: the original aa sequence of a peptide (base_seq),
            and the relative position of the modification.
    output: returns the list of all codons possibly associated 
            with the substitution, presented as a string separated
            by white spaces.
    """
    possible_codons = []
    proteins = find_proteins(base_seq, names_list, boundaries_aa, W_aa, sa)
    if proteins:
        proteins = proteins.split(" ")
    else:
        return '_'
    for p in proteins:
        if p in record_dict:
            s = record_dict[p].seq
            seq_i = s.translate().find(base_seq)
            i = seq_i + modified_pos
            possible_codons.append(codonify(s)[i])
        else:
            possible_codons.append('_')
    return " ".join(possible_codons)


def fetch_best_codons(modified_seq, record_dict, names_list, boundaries_aa, W_aa, sa):
    """
    input: a modified sequence, e.g. LQV(0.91)A(0.09)EK
    output: the list of codons associated with the most likely
            position
    """
    possible_sites = re.findall('\(([^\)]+)\)', modified_seq)
    best_site = np.argmax([float(i) for i in possible_sites])
    modified_pos_prime = [m.start() - 1 for m in re.finditer('\(', modified_seq)][best_site]
    modified_pos = len(re.sub('\(([^\)]+)\)', '', modified_seq[:modified_pos_prime]))
    base_seq = re.sub('\(([^\)]+)\)', '', modified_seq)
    return fetch_codon(base_seq, modified_pos, record_dict, names_list, boundaries_aa, W_aa, sa)


def find_substitution_position_local(modified_seq, protein, record_dict):
    """
    returns the position of a substitutions relative to the start
    of the protein sequence
    """
    possible_sites = re.findall('\(([^\)]+)\)', modified_seq)
    best_site = np.argmax([float(i) for i in possible_sites])
    modified_pos_prime = [m.start() - 1 for m in re.finditer('\(', modified_seq)][best_site]
    modified_pos = len(re.sub('\(([^\)]+)\)', '', modified_seq[:modified_pos_prime]))
    base_seq = re.sub('\(([^\)]+)\)', '', modified_seq)
    s = record_dict[protein].seq
    seq_i = s.translate().find(base_seq)
    i = seq_i + modified_pos
    return i


def find_positions_local(modified_seq, proteins, record_dict):
    """
    returns the position of a substitutions relative to the start
    of the protein sequence, across all the codons
    """
    positions = []
    for prot in proteins.split(" "):
        positions.append(str(find_substitution_position_local(modified_seq, prot, record_dict)))
    return " ".join(positions)


def protein_lengths(proteins, record_dict):
    """
    @author: Tehila
    @return: the length of proteins with the substitutions
    """
    lengths = []
    for prot in proteins.split(" "):
        s = record_dict[prot].seq
        lengths.append(str(len(s)))
    return " ".join(lengths)


def is_gene(record):
    """annotation by Tehila: validate CDS by making sure it's dividing by 3,
    starting with start codon, including non-canonical start codons, and end with end codon"""
    if len(record.seq) % 3 != 0:
        return False
    if not record.seq[:3] in {'ATG', 'GTG', 'TTG', 'ATT', 'CTG'}:  # TEHILA: why 5 different options???
        return False
    if record.seq[-3:].translate() != '*':
        return False
    return True


def count_non_uppercase(string):
    non_uppercase = re.findall(r'[^A-Z]', string)
    return len(non_uppercase)


def refine_localization_probabilities(modified_seq, threshold=AA_LOC_THRESH):
    """
    returns the AAs that were possibly modified (with p > threshold).
    Input: modified sequence (a string of AA with p of each to contain modification: APKV(.7)ML(.3)L means that V was modified with p = .7 and L with p = .3)
    Output: A string with all candidate AAs separated by ';' (V;L).
    """
    modified_sites = [(m.start() - 1, modified_seq[m.start() - 1]) for m in re.finditer('\(', modified_seq)]
    weights = [float(i) for i in re.findall('\(([^\)]+)\)', modified_seq)]
    if max(weights) > threshold:
        position, aa = modified_sites[np.argmax(weights)]
        position = position - count_non_uppercase(modified_seq[:position])
        return position, aa
    else:
        return None, None


def peptide_termini(peptide_len, error_position):
    if error_position - 1 == peptide_len:
        return 'c-term'
    if error_position == 0:
        return 'n-term'


def protein_termini(protein_len, error_position):
    if error_position - 1 == protein_len:
        return 'c-term'
    if error_position == 0:
        return 'n-term'


def is_prot_nterm(sequence, W_aa, sa):
    """
    Does the peptide originate at the protein's N-term
    """
    for start in SA_search(sequence, W_aa, sa):
        if W_aa[start - 1] == '*':
            return True
        if W_aa[start - 2] == '*':
            return True
    return False


def get_mispairing_mask():
    mask = pd.DataFrame(data=False, index=CODONS, columns=list('ACDEFGHKLMNPQRSTVWY'), dtype=float)
    for label in CODONS:
        near_cognates = np.array([hamming(i, label) == 1 for i in CODONS])
        reachable_aa = set(np.array(list(AMINO_ACIDS))[near_cognates])
        mask.loc[label] = [i in reachable_aa for i in 'ACDEFGHKLMNPQRSTVWY']

    inverted_codon_table = get_inverted_codon_table()
    for label in mask.index:  # removes "near-cognates" that encodes the same AA
        for col in mask.columns:
            if label in inverted_codon_table[col]:
                mask.loc[label, col] = float('NaN')
    return mask


def hamming(s1, s2): return sum(a != b for a, b in zip(s1, s2))


def is_mispairing(row, mask):
    """
    Returns whether the substitution is mispairing or misloading, based on the
    near-cognate mask.
    """
    codon = row['codon']
    destination = row['destination']
    if pd.notnull(codon) and pd.notnull(destination):
        if (codon in mask.index) and destination:
            return mask.loc[codon, destination]
        else:
            return 0
    else:
        return float('NaN')


def suffix_array(text, _step=16):
    # TEHILA: suffix array is a data structure that help to find
    # substrings or patterns in text. LCP is an auxiliary data structure to the suffix array.
    # look for it on Wikipedia
    """Analyze all common strings in the text.
    
    Short substrings of the length _step a are first pre-sorted. The are the 
    results repeatedly merged so that the garanteed number of compared
    characters bytes is doubled in every iteration until all substrings are
    sorted exactly.
    
    Arguments:
        text:  The text to be analyzed.
        _step: Is only for optimization and testing. It is the optimal length
               of substrings used for initial pre-sorting. The bigger value is
               faster if there is enough memory. Memory requirements are
               approximately (estimate for 32 bit Python 3.3):
                   len(text) * (29 + (_size + 20 if _size > 2 else 0)) + 1MB
    
    Return value:      (tuple)
      (sa, rsa, lcp)
        sa:  Suffix array                  for i in range(1, size):
               assert text[sa[i-1]:] < text[sa[i]:]
        rsa: Reverse suffix array          for i in range(size):
               assert rsa[sa[i]] == i
        lcp: Longest common prefix         for i in range(1, size):
               assert text[sa[i-1]:sa[i-1]+lcp[i]] == text[sa[i]:sa[i]+lcp[i]]
               if sa[i-1] + lcp[i] < len(text):
                   assert text[sa[i-1] + lcp[i]] < text[sa[i] + lcp[i]]
    >>> suffix_array(text='banana')
    ([5, 3, 1, 0, 4, 2], [3, 2, 5, 1, 4, 0], [0, 1, 3, 0, 0, 2])
    
    Explanation: 'a' < 'ana' < 'anana' < 'banana' < 'na' < 'nana'
    The Longest Common String is 'ana': lcp[2] == 3 == len('ana')
    It is between  tx[sa[1]:] == 'ana' < 'anana' == tx[sa[2]:]
    """
    tx = text
    size = len(tx)
    step = min(max(_step, 1), len(tx))
    sa = list(range(len(tx)))
    sa.sort(key=lambda i: tx[i:i + step])
    grpstart = size * [False] + [True]  # a boolean map for iteration speedup.
    # It helps to skip yet resolved values. The last value True is a sentinel.
    rsa = size * [None]
    stgrp, igrp = '', 0
    for i, pos in enumerate(sa):
        st = tx[pos:pos + step]
        if st != stgrp:
            grpstart[igrp] = (igrp < i - 1)
            stgrp = st
            igrp = i
        rsa[pos] = igrp
        sa[i] = pos
    grpstart[igrp] = (igrp < size - 1 or size == 0)
    while grpstart.index(True) < size:
        # assert step <= size
        nextgr = grpstart.index(True)
        while nextgr < size:
            igrp = nextgr
            nextgr = grpstart.index(True, igrp + 1)
            glist = []
            for ig in range(igrp, nextgr):
                pos = sa[ig]
                if rsa[pos] != igrp:
                    break
                newgr = rsa[pos + step] if pos + step < size else -1
                glist.append((newgr, pos))
            glist.sort()
            for ig, g in groupby(glist, key=itemgetter(0)):
                g = [x[1] for x in g]
                sa[igrp:igrp + len(g)] = g
                grpstart[igrp] = (len(g) > 1)
                for pos in g:
                    rsa[pos] = igrp
                igrp += len(g)
        step *= 2
    del grpstart
    del rsa
    return sa


def SA_search(P, W, sa):
    """
    W - a string (of concatenated peptides)
    sa - a suffix array for W
    P - pattern to look for in sa
    """
    lp = len(P)
    n = len(sa)
    l = 0;
    r = n
    while l < r:
        mid = (l + r) // 2
        a = sa[mid]
        if P > W[a: a + lp]:
            l = mid + 1
        else:
            r = mid
    s = l;
    r = n
    while l < r:
        mid = (l + r) // 2
        a = sa[mid]
        if P < W[a: a + lp]:
            r = mid
        else:
            l = mid + 1
    return [sa[i] for i in range(s, r)]


# TEHILA: my improvement:
def find_homologous_peptide(P, W_aa_ambiguous, sa_ambiguous, trypsin_digested):
    """
    Gets a peptide and returns whether it has homolegous in the genome.
    If so, that peptide is discarded.
    It have two modes: if the peptides are trypsin digested, it will only
    look for peptides that could have been made by trypsin digestion.
    otherwise it would assume that the peptide could be everywhere
    """
    if trypsin_digested:
        if len(SA_search('K' + P, W_aa_ambiguous, sa_ambiguous)) > 0:
            return False
        if len(SA_search('R' + P, W_aa_ambiguous, sa_ambiguous)) > 0:
            return False
        if len(SA_search('*' + P, W_aa_ambiguous, sa_ambiguous)) > 0:
            return False
        return True
    else:
        if len(SA_search(P, W_aa_ambiguous, sa_ambiguous)) > 0:
            return False
        return True


# def create_modified_seq(modified_seq, destination):
#     possible_sites = re.findall('\(([^\)]+)\)', modified_seq)
#     best_site = np.argmax([float(i) for i in possible_sites])
#     modified_pos_prime = [m.start() - 1 for m in re.finditer('\(', modified_seq)][best_site]
#     modified_pos = len(re.sub('\(([^\)]+)\)', '', modified_seq[: modified_pos_prime]))
#     base_seq = re.sub('\(([^\)]+)\)', '', modified_seq)
#     return base_seq[: modified_pos] + destination + base_seq[modified_pos + 1:]

def create_modified_seq(base_seq: str, err_pos_pep: float, destination: str, xle='L'):
    base_seq_list = list(base_seq)
    if destination == 'L':
        destination = xle
    base_seq_list[int(err_pos_pep)] = destination
    return ''.join(base_seq_list)

def ask_about_digestion():
    answer = input("is the sample trypsin-digested?\nIf yes, press enter.\nOtherwise, for not known digestion,"
                   "press 1 and then enter.\nFor further information check the function find_homologous_peptide"
                   " at subsUtils.py file")
    if answer == '':
        print("Keep processing as trypsin-digested data")
        return True
    elif answer == "1":
        print("Keep processing as unspecific-digested data")
        return False
    else:
        raise SystemExit('Your input is not valid')


def sum_intensities(intensities):
    inten_list = intensities.split(';')
    return sum(map(float, inten_list))


def fasta2string(fasta_path=PEP_FASTA):
    records = list(SeqIO.parse(fasta_path, "fasta"))
    fasta_str = ''.join(map(lambda rec: str(rec.seq)+'*', records))
    return fasta_str