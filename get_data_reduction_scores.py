#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 07:58:44 2024

@author: charmainechia
"""
import os
import numpy as np
from utils.utils import fetch_sequences_from_fasta, write_sequence_to_fasta, get_ref_seq_idxs_aa_from_msa, compute_entropy

def ShanEntropy(msa_fname, output_fname, ref_seq_idxs, ref_seq):
    from data_reduction.alfa2cons import alfa2cons
    csv = alfa2cons('./data/msa/'+msa_fname+'.fasta', './data/data_reduction/'+output_fname+'.csv', save_csv=False)
    csv_filt = csv[csv['Position'].isin(ref_seq_idxs)].copy()
    ref_seq_parsed = list(ref_seq) if ref_seq is not None else None
    if ref_seq_parsed is not None:
        csv_filt['AA'] = ref_seq_parsed
    csv_filt['RealPos'] = list(np.arange(len(ref_seq_parsed)) + 1)
    csv_filt = csv_filt.rename(columns={'Position': 'PositionMSA'})
    csv_filt = csv_filt[['AA', 'RealPos', 'PositionMSA', 'shanID', 'shanMS']]
    # write dataframe to file
    csv_filt.to_csv('./data/data_reduction/'+output_fname+'_ShanEntropy.csv')

def SIFT(msa_fname, output_fname, ref_seq_name):
    from data_reduction import access_sift_webserver
    msa_path = './data/msa/'+msa_fname+'.fasta'
    msa_seqs, msa_names, _ = fetch_sequences_from_fasta(msa_path)

    msa_idx = msa_names.index(ref_seq_name)
    msa_names_rearranged = [msa_names[msa_idx]] + msa_names[:msa_idx] + msa_names[msa_idx+1:]
    msa_seqs_rearranged = [msa_seqs[msa_idx]] + msa_seqs[:msa_idx] + msa_seqs[msa_idx+1:]
    msa_path_rearranged = './data/msa/'+msa_fname + f'_{ref_seq_name}.fasta'
    write_sequence_to_fasta(msa_seqs_rearranged, msa_names_rearranged, msa_fname+f'_{ref_seq_name}', './data/msa/')
    csv = access_sift_webserver.main(
      {'-a': os.path.abspath(msa_path_rearranged),
       '--fname': './data/data_reduction/'+msa_fname+f'_sift_{ref_seq_name}.csv'
       }
    )
    # get probabilities pre-normalization
    probs_matrix_norm = np.transpose(csv.iloc[:, -20:].to_numpy()).astype(float)
    csv.iloc[:,-20:] = np.transpose(probs_matrix_norm)
    prob_col = csv['prob'].to_numpy()
    probs_matrix = probs_matrix_norm * prob_col
    # get average sift score for each residue
    probs_mean = np.mean(probs_matrix_norm, axis=0)
    csv.insert(2, 'sift_norm', list(probs_mean))
    # calculate entropy for each position
    ent = compute_entropy(probs_matrix)
    csv.insert(2, 'entropy', list(ent))
    # write dataframe to file
    csv = csv.rename(columns={'pos': 'RealPos', 'wt': 'AA'})
    csv.to_csv('./data/data_reduction/'+output_fname+'_sift.csv')
    os.remove(msa_path_rearranged)
    os.remove('./data/data_reduction/'+msa_fname+f'_sift_{ref_seq_name}.csv')


def main():
    msa_fname = 'GOh1052_msa'
    ref_seq_name = 'FGGALOX'
    ref_seq = None
    output_fname = 'GOh1052'
    if ref_seq is None:
        _, ref_seq, ref_seq_idxs = get_ref_seq_idxs_aa_from_msa('./data/msa/'+msa_fname+'.fasta', ref_seq_name)

    # run ShanEntropy
    ShanEntropy(msa_fname, output_fname, ref_seq_idxs, ref_seq_name, ref_seq)
    # run SIFT
    SIFT(msa_fname, output_fname, ref_seq_name)

if __name__ == "__main__":
    main()