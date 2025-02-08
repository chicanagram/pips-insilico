#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 07:58:44 2024

@author: charmainechia
"""
try:
    import pandas as pd
    pandas_imported = True
except ImportError as e:
    pandas_imported = False
import numpy as np
import os
import platform
opsys = platform.system()

def sort_list(lst):
    lst.sort()
    return lst

def intersection(lst1, lst2):
    return list(set(lst1) & set(lst2))

def get_mutstr(mutation):
    if isinstance(mutation, list):
        mutstr = '+'.join(mutation)
    else:
        mutstr = mutation
        mutation = [mutation]
    return mutstr, mutation

def get_mapping(mapping, res):
    # Map residue initial to code
    result = mapping[res]
    return result

def mkDir(res, output_dir, remove_existing_dir=True):
    import shutil
    # making new directory
    new_dir = (output_dir + res)
    if os.path.exists(new_dir):
        # remove if directory exists, and make new directory
        if remove_existing_dir:
            shutil.rmtree(new_dir)
            os.makedirs(new_dir)
    else:
        os.makedirs(new_dir)
    return new_dir

def findProcess(process_name):
    if opsys=='Windows':
        return [int(item.split()[1]) for item in os.popen('tasklist').read().splitlines()[4:] if process_name in item.split()]
    elif opsys=='Linux' or opsys=='Darwin':
        return [int(pid) for pid in os.popen('pidof '+process_name).read().strip(' \n').split(' ')]

def exit_program(pid):
    import signal
    print("Sending SIGINT to self...")
    os.kill(pid, signal.SIGINT)
    print('Exited program', pid)
def is_float(string):
    if string.replace(".", "").replace("-", "").isnumeric():
        return True
    else:
        return False

def compute_entropy(probs_matrix):
    """
    Computes the entropy for each position in the given probability matrix.
    """
    from scipy.stats import entropy

    # Initialize an empty list to store the entropy values
    entropy_values = []

    # Iterate over the columns of probs_matrix
    for i in range(probs_matrix.shape[1]):
        # Compute the entropy for the probabilities at the current position
        H = entropy(probs_matrix[:, i], base=2)
        entropy_values.append(H)

    # Convert entropy_values to a numpy array for convenience
    entropy_values = np.array(entropy_values)
    return entropy_values

def save_dict_as_csv(datadict, cols, log_fpath, csv_suffix ='', multiprocessing_proc_num=None):
    # save results as CSV
    csv_txt = ''
    # get csv_suffix if running multiprocessing
    if multiprocessing_proc_num is not None:
        csv_suffix += '_' + str(multiprocessing_proc_num)

    # check if file exists yet
    log_fpath_full = log_fpath + csv_suffix + '.csv'
    if not os.path.exists(log_fpath_full):
        # if not, start a new file with headers
        write_mode = 'w'
        csv_txt += ','.join(cols) + '\n'
    else:
        write_mode = 'a'

    # convert dict of lists to list of dicts
    if isinstance(datadict[cols[0]], list):
        num_rows = len(datadict[cols[0]])
        datadict_byrow = []
        for row_idx in range(num_rows):
            row = []
            for col in cols:
                row.append(datadict[col][row_idx])
            datadict_byrow.append(row)
    else:
        row = []
        for col in cols:
            row.append(datadict[col])
        datadict_byrow = [row]

    # add data to csv file
    for row in datadict_byrow:
        csv_txt += ','.join([str(el) for el in row])
        csv_txt += '\n'
    # save the changes
    with open(log_fpath_full, write_mode) as f:
        f.write(csv_txt)
    return csv_txt, log_fpath_full, write_mode

def combine_csv_files(log_fpath_list, output_dir, output_fname, remove_combined_files=True):
    # combine files spawned
    txt_all_list = []
    for i, log_fpath in enumerate(log_fpath_list):
        with open(log_fpath, 'r') as f:
            if i==0:
                txt_all_list += f.readlines()
            else:
                txt_all_list += f.readlines()[1:]
    if os.path.exists(output_dir + output_fname + '.csv'):
        write_mode = 'a'
        txt_all_list = txt_all_list[1:]
    else:
        write_mode = 'w'
    # get text string to write
    txt_all = '\n'.join(txt_all_list)
    txt_all = txt_all.replace('\n\n', '\n').replace(',\n', '\n')
    # update or save file
    with open(output_dir + output_fname + '.csv', write_mode) as f:
        f.write(txt_all)
    # remove combined files
    if remove_combined_files:
        for log_fpath in log_fpath_list:
            os.remove(log_fpath)
    return txt_all

def split_mutation(mutation, aa_letter_representation=False):
    # Convert point mutation to wildtype residue, muted residue and mutation position
    mutation = list(mutation)
    WT_res = mutation[0]
    MUT_res = mutation[-1]
    if not aa_letter_representation:
        WT_res = get_mapping(mapping, WT_res)
        MUT_res = get_mapping(mapping, MUT_res)
    MUT_pos = mutation[1:len(mutation)-1]
    MUT_pos = int(''.join(MUT_pos))
    return WT_res, MUT_pos, MUT_res
def split_wildtype(mutation):
    # Convert point mutation to wildtype residue, muted residue and mutation position
    WT_res = mutation[0]
    MUT_pos = int(mutation[1:len(mutation)-1])
    MT_res = mutation[-1]
    return WT_res, MUT_pos, MT_res
def get_mutations(wildtype_list):
    # get amino acid list to perform mutations to
    aaList = ['A', 'H', 'Y', 'R', 'T', 'K', 'M', 'D', 'N', 'C', 'Q', 'E', 'G', 'I', 'L', 'F', 'P', 'S', 'W', 'V']
    # get all mutations to run
    mutations = []
    for wt in wildtype_list:
        wtAA = wt[0]
        for aa in aaList:
            if aa != wtAA:
                mt = wt + aa
                mutations.append(mt)
    print('mutant:', mutations)
    return mutations

def get_mutated_sequence(seq_base, mutations, seq_name=None, write_to_fasta=None):
    mutations_list = []
    seq_name_list = []
    sequence_list = []
    fasta_list = []

    # get mutated sequences
    for mut in mutations:
        if mut is not None:
            wildtype_aa, position, mutant_aa = split_wildtype(mut)

            if seq_base[position-1] == wildtype_aa:
                if seq_name is None:
                    seq_name_wmut = mut
                else:
                    seq_name_wmut = seq_name + '_' + mut
            else:
                seq_name_wmut = seq_name + '_' + mut
            list_seq_base = list(seq_base)
            list_seq_base[position-1] = mutant_aa
            seq_mutate = "".join(list_seq_base)

        else:
            mut = 'WT'
            seq_mutate = seq_base
            seq_name_wmut = seq_name + '_' + mut

        mutations_list.append(mut)
        sequence_list.append(seq_mutate)
        seq_name_list.append(seq_name_wmut)

        # write to fasta
        if write_to_fasta is not None:
            fasta_file = write_sequence_to_fasta(seq_mutate, seq_name_wmut, seq_name_wmut, write_to_fasta)
            fasta_list.append(fasta_file)

    return mutations_list, seq_name_list, sequence_list, fasta_list
def get_enz_shortname(enzname):
    if enzname in ['CviUPO', 'DcaUPO', 'TruUPO', 'HspUPO']:
        return enzname
    else:
        return enzname[:2]+enzname[-3:]

def get_enzname_from_expdata(csv_fpath, df=None):
    if df is None:
        df = pd.read_csv(csv_fpath, index_col=0)
    struct_list = df.struct.tolist()
    enz_list = ['_'.join(struct.split('_')[:-1]) for struct in struct_list]
    enz_short_list = [get_enz_shortname(enz) for enz in enz_list]
    df['enz_name'] = enz_list
    df['enz_name_short'] = enz_short_list
    print(df)
    df.to_csv(csv_fpath)
    unique_enz_list = list(set(enz_list))
    unique_enzshort_list = list(set(enz_short_list))
    return df, unique_enz_list, unique_enzshort_list
def read_fasta(sequence_fpath):
    seq_dict = {}
    with open(sequence_fpath) as f:
        lines = f.readlines()
    for i in range(0, len(lines)):
        s = lines[i]
        if s[0] == '>':
            seq_name = s[1:].strip()
        else:
            seq_dict[seq_name] = s.strip()
    return seq_dict
def fetch_sequences_from_fasta(sequence_fpath):
    from Bio import SeqIO
    import os
    sequence_names = []
    sequence_list = []
    sequence_descriptions = []
    for j, record in enumerate(SeqIO.parse(sequence_fpath, "fasta")):
        sequence_names.append(record.id)
        sequence_list.append(str(record.seq))
        sequence_descriptions.append(record.description)
    return sequence_list, sequence_names, sequence_descriptions
def write_sequence_to_fasta(sequences, seq_names, filename, fasta_dir):
    fasta_file = fasta_dir + filename + '.fasta'
    if isinstance(sequences, str) and isinstance(seq_names, str):
        sequences = [sequences]
        seq_names = [seq_names]
    with open(fasta_file, 'w') as f:
        for i, (sequence, seq_name) in enumerate(zip(sequences, seq_names)):
            f.write('> ' + seq_name + '\n')
            if i==len(sequences)-1:
                f.write(sequence)
            else:
                f.write(sequence+'\n')
    print('Saved fasta file to ' + fasta_file)
    return fasta_file

def split_fasta(max_res_per_fasta, fasta_fname, fasta_dir):
    seqs, seq_names, _ = fetch_sequences_from_fasta(fasta_dir+fasta_fname)
    num_seq = len(seqs)
    seq_len = len(seqs[0])
    print('# of seqs:', num_seq, '; seq len:', seq_len)
    num_seq_per_fasta = int(np.floor(max_res_per_fasta / seq_len))
    num_fasta_to_split = int(np.floor(num_seq / num_seq_per_fasta))
    print('# fasta files to split too:', num_fasta_to_split, '; # of seqs per fasta:', num_seq_per_fasta)

    for i in range(num_fasta_to_split):
        seqs_i = seqs[i * num_seq_per_fasta:min((i + 1) * num_seq_per_fasta, num_seq)]
        seq_names_i = seq_names[i * num_seq_per_fasta:min((i + 1) * num_seq_per_fasta, num_seq)]
        fasta_out = write_sequence_to_fasta(seqs_i, seq_names_i, fasta_fname.replace('.fasta','') + '_' + str(i), fasta_dir)

def deduplicate_fasta_sequences(fasta_fname_in, fasta_fname_out=None, fasta_dir='./'):
    seqs, seq_names, _ = fetch_sequences_from_fasta(fasta_dir+fasta_fname_in)
    seqs_deduped = []
    seq_names_deduped = []
    for seq, seq_name in zip(seqs, seq_names):
        if seq_name not in seq_names_deduped:
            seq_names_deduped.append(seq_name)
            seqs_deduped.append(seq)
    print('# of seqs [BEFORE]', len(seq_names), '[AFTER]', len(seq_names_deduped))
    # save resulting fasta
    if fasta_fname_out is None:
        fasta_fname_out = fasta_fname_in
    fasta_fname_out = fasta_fname_out.replace('.fasta','')
    fasta_out = write_sequence_to_fasta(seqs_deduped, seq_names_deduped, fasta_fname_out,  fasta_dir)
    return fasta_out, seqs_deduped, seq_names_deduped

def get_ref_seq_idxs_aa_from_msa(msa_path, ref_seq_names, zero_indexed=False):
    if isinstance(ref_seq_names, str):
        ref_seq_name_list = [ref_seq_names]
    else:
        ref_seq_name_list = ref_seq_names

    # get ref_seq_name_list found in MSA
    ref_seq_name_list_inmsa = []
    ref_seq_idxs_list_inmsa = []
    ref_seq_list_inmsa = []
    # check that all ref seqs are in the MSA
    msa_seqs, msa_names, _ = fetch_sequences_from_fasta(msa_path)
    for ref_seq_name in ref_seq_name_list:
        ref_seq_inmsa = ref_seq_name in msa_names
        print(ref_seq_name + ' is in MSA: ' + str(ref_seq_inmsa))
        if ref_seq_name in msa_names:
            ref_seq_name_list_inmsa.append(ref_seq_name)
            msa_idx = msa_names.index(ref_seq_name)
            msa_seq = msa_seqs[msa_idx]
            seq_filt = ''
            idx_filt = []
            for i, letter in enumerate(list(msa_seq)):
                if letter != '-':
                    seq_filt += letter
                    if zero_indexed:
                        idx_filt.append(i)
                    else:
                        idx_filt.append(i+1)
            ref_seq_idxs_list_inmsa.append(idx_filt)
            ref_seq_list_inmsa.append(seq_filt)

    if isinstance(ref_seq_names, str):
        ref_seq_name_list_inmsa = ref_seq_name_list_inmsa[0]
        ref_seq_list_inmsa = ref_seq_list_inmsa[0]
        ref_seq_idxs_list_inmsa = ref_seq_idxs_list_inmsa[0]

    return ref_seq_name_list_inmsa, ref_seq_list_inmsa, ref_seq_idxs_list_inmsa

def add_sequences_to_dataframe(csv_fpath, enzname_to_seq_dict, df=None):
    # get dataframe
    if df is None:
        df = pd.read_csv(csv_fpath, index_col=0)
    # get enzyme names
    enz_name_list = df['enz_name'].tolist()
    # match to sequence
    sequence_list = []
    for enz_name in enz_name_list:
        sequence_list.append(enzname_to_seq_dict[enz_name])
    # update df
    df['sequence'] =  sequence_list
    # save df
    df.to_csv(csv_fpath)
    return df
def get_sequences_for_enzyme_list(csv_fpath, sequence_fpath, create_new_fasta=None, df=None):
    # get dataframe and enzyme lists
    df, unique_enz_list, unique_enzshort_list = get_enzname_from_expdata(csv_fpath, df=df)

    # get sequences
    seq_list, seq_names, sequence_descriptions = fetch_sequences_from_fasta(sequence_fpath)
    print(sequence_descriptions)

    # match enzymes in list from dataframe to sequences in sequence file
    enzname_to_seq_dict = {}
    for enz in unique_enz_list:
        if enz in seq_names:
            match_idx = seq_names.index(enz)
            enzname_to_seq_dict.update({enz:seq_list[match_idx]})
            continue
        elif enz not in seq_names:
            try:
                match_idx = [i for i, des in enumerate(sequence_descriptions) if enz in des][0]
                enzname_to_seq_dict.update({enz:seq_list[match_idx]})
            except:
                print('No match found for ' + enz + ', ')

    unmatched_enz_list = [enz for enz in unique_enz_list if enz not in enzname_to_seq_dict]
    print('Matched enzymes:', len(enzname_to_seq_dict), enzname_to_seq_dict)
    print('Unmatched enzymes:', len(unmatched_enz_list), unmatched_enz_list)

    # save enzymes and sequences to fasta
    if create_new_fasta is not None:
        write_sequence_to_fasta(list(enzname_to_seq_dict.values()), list(enzname_to_seq_dict.keys()), create_new_fasta, os.path.dirname(sequence_fpath)+'/')

    # update dataframe with sequences
    df = add_sequences_to_dataframe(csv_fpath, enzname_to_seq_dict)

    return enzname_to_seq_dict, unique_enz_list, unique_enzshort_list, unmatched_enz_list, df

def extract_yasara_logdata_to_csv(data_folder, dataset_label, mutations_present=False, res_type_list=['fullbreakdown']):
    yasara_output_folder = data_folder + yasara_subfolder + 'Output/' + dataset_label + '/'
    DDG_cols_nomut = {
        'fullbreakdown': [
            'avg_epotWTRtr', 'avg_esolWTcolRtr', 'avg_esolWTvdwRtr', 'avg_surfaccWTRtr',
            'avg_epotWTCpx', 'avg_esolWTcolCpx', 'avg_esolWTvdwCpx', 'avg_surfaccWTCpx',
            'avg_epotWTLgd', 'avg_esolWTcolLgd', 'avg_esolWTvdwLgd', 'avg_surfaccWTLgd',
            'avg_surfcost', 'avg_ebindWT', 'ligname', 'setname', 'struct'
        ],
        'avg': ['struct', 'setname', 'ligname', 'avg_ebindWT']
    }
    DDG_cols_mut = {
        'fullbreakdown': [
            'avg_epotWTRtr', 'avg_esolWTcolRtr', 'avg_esolWTvdwRtr', 'avg_surfaccWTRtr',
            'avg_epotMTRtr', 'avg_esolMTcolRtr', 'avg_esolMTvdwRtr', 'avg_surfaccMTRtr',
            'avg_epotWTCpx', 'avg_esolWTcolCpx', 'avg_esolWTvdwCpx', 'avg_surfaccWTCpx',
            'avg_epotMTCpx', 'avg_esolMTcolCpx', 'avg_esolMTvdwCpx', 'avg_surfaccMTCpx',
            'avg_epotWTLgd', 'avg_esolWTcolLgd', 'avg_esolWTvdwLgd', 'avg_surfaccWTLgd',
            'avg_surfcost', 'avg_DDG', 'mutname', 'ligname', 'setname', 'struct'
        ],
        'avg': ['struct', 'setname', 'ligname', 'avg_DDG']
    }
    logfiles = {
        'fullbreakdown': [f for f in os.listdir(yasara_output_folder) if (f.endswith('log') and f.startswith('DDG_fullbreakdown'))],
        'avg': [f for f in os.listdir(yasara_output_folder) if (f.endswith('log') and f.startswith('DDG_avg'))]
    }
    if mutations_present:
        DDG_cols = DDG_cols_mut.copy()
    else:
        DDG_cols = DDG_cols_nomut.copy()
    binding_res = {}
    for res_type in res_type_list:
        binding_res[res_type] = []
        logfiles_byrestype = logfiles[res_type]
        print(res_type)
        for logfile in logfiles_byrestype:
            with open(yasara_output_folder+logfile, 'r') as f:
                lines = f.readlines()
                for line in lines:
                    line = line.rstrip()
                    vals = line.split(' ')
                    vals = [float(val) if is_float(val) else val for val in vals ]
                    d = {}
                    for k, v in zip(DDG_cols[res_type], vals): d.update({k:v})
                    binding_res[res_type].append(d)
        binding_res[res_type] = pd.DataFrame(binding_res[res_type])[DDG_cols[res_type]]
        binding_res[res_type].to_csv(yasara_output_folder + 'DDG_' + res_type + '_' + dataset_label + '.csv')
        print(binding_res[res_type].shape)
        print(binding_res[res_type])
    return binding_res

######################
## POTTS MODEL UTILS ##
######################

def parse_fasta_to_numpy(filename):
    '''function to parse fasta file'''
    header = []
    sequence = []
    lines = open(filename, "r")
    for line in lines:
        line = line.rstrip()
        if line[0] == ">":
            header.append(line[1:])
            sequence.append([])
        else:
            sequence[-1].append(line)
    lines.close()
    sequence = [''.join(seq) for seq in sequence]
    return np.array(header), np.array(sequence)

def one_hot(msa, states):
    one = np.eye(states)
    return one[msa]

def mk_msa(seqs):
    '''one hot encode msa'''
    ################
    alphabet = "ARNDCQEGHILKMFPSTWYV-"
    states = len(alphabet)
    a2n = {}
    for a, n in zip(alphabet, range(states)):
        a2n[a] = n
    def aa2num(aa):
        '''convert aa into num'''
        if aa in a2n:
            return a2n[aa]
        else:
            return a2n['-']
    ################
    msa = []
    for seq in seqs:
        msa.append([aa2num(aa) for aa in seq])
    msa_ori = np.array(msa)
    return msa_ori, one_hot(msa_ori, states)

def get_mtx(W):
    # l2norm of 20x20 matrices (note: we ignore gaps)
    raw = np.sqrt(np.sum(np.square(W[:, :, :, :]), (1, 3)))
    np.fill_diagonal(raw, 0)
    # apc (average product correction)
    ap = np.sum(raw, 0, keepdims=True) * np.sum(raw, 1, keepdims=True) / np.sum(raw)
    apc = raw - ap
    np.fill_diagonal(apc, 0)
    return (raw, apc)
    ################
    msa = []
    for seq in seqs:
        msa.append([aa2num(aa) for aa in seq])
    msa_ori = np.array(msa)
    return msa_ori, one_hot(msa_ori, states)

def get_contacts(model, save_res=None):
    # get trained weights
    w = model.GREMLIN_.W0.detach().numpy()
    w = (w + w.T).reshape(L, A, L, A)
    raw, apc = get_mtx(w)
    # get contacts
    if save_res is not None:
        np.savetxt(save_res + '_raw.csv', raw, delimiter=",")
        np.savetxt(save_res + '_apc.csv', apc, delimiter=",")
    return raw, apc

def annotate_heatmap(array_2D, ax, ndecimals=2, fontsize=8):
    for (j,i),label in np.ndenumerate(array_2D):
        if ndecimals==0 and ~np.isnan(label):
            label = int(label)
        else:
            label = round(label,ndecimals)
        ax.text(i,j,label,ha='center',va='center', color='0.8', fontsize=fontsize, fontweight='bold')

def symlog(data):
    idx_pos = np.where(data > 0)
    idx_neg = np.where(data < 0)
    data_symlog = np.zeros((data.shape[0], data.shape[1]))
    data_symlog[:] = np.nan
    data_symlog[idx_pos] = np.log(data[idx_pos])
    data_symlog[idx_neg] = -np.log(-data[idx_neg])
    return data_symlog

def heatmap(array, c='viridis', ax=None, cbar_kw={}, cbarlabel="", datamin=None, datamax=None, logscale_cmap=False,
            annotate=None, row_labels=None, col_labels=None, show_gridlines=True, fontsize=8):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (N, M).
    row_labels
        A list or array of length N with the labels for the rows.
    col_labels
        A list or array of length M with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    """
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    if not ax:
        ax = plt.gca()
    cmap = getattr(plt.cm, c)

    # get array size and xy labels
    data = array.astype(float)
    ny, nx = data.shape

    # get row and column labels
    if row_labels is None:
        row_labels = list(np.arange(ny) + 1)
    if col_labels is None:
        col_labels = list(np.arange(nx) + 1)

    # get locations of nan values and negative values, replace values so these don't trigger an error
    naninds = np.where(np.isnan(data) == True)
    infinds = np.where(np.isinf(data) == True)
    if len(infinds[0]) > 0:
        data[infinds] = np.nan
    if len(naninds[0]) > 0:
        data[naninds] = np.nanmean(data)
    if len(infinds[0]) > 0:
        data[infinds] = np.nanmean(data)
    data_cmap = data.copy()

    # get min and max values
    if datamin is None:
        datamin = np.nanmin(data_cmap)
    if datamax is None:
        datamax = np.nanmax(data_cmap)

    # get colormap to plot
    if logscale_cmap:  # plot on logscale
        data_cmap = symlog(data_cmap)
        datamin, datamax = np.min(data_cmap), np.max(data_cmap)

    # get cmap gradations
    dataint = (datamax - datamin) / 100
    norm = plt.Normalize(datamin, datamax + dataint)
    # convert data array into colormap
    colormap = cmap(norm(data_cmap))

    # Set the positions of nan values in colormap to 'lime'
    colormap[naninds[0], naninds[1], :3] = 0, 1, 0
    colormap[infinds[0], infinds[1], :3] = 1, 1, 1

    # plot colormap
    im = ax.imshow(colormap, interpolation='nearest')

    # Create colorbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.07)
    cbar = ax.figure.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax)

    if logscale_cmap == True:
        cbar_labels = cbar.ax.get_yticks()
        cbar.set_ticks(cbar_labels)
        cbar_labels_unlog = list(np.round(np.exp(np.array(cbar_labels)), 2))
        cbar.set_ticklabels(cbar_labels_unlog)

    # Turn off gridlines if required
    ax.tick_params(axis='both', which='both', length=0, gridOn=show_gridlines)

    # We want to show all ticks...
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    # ... and label them with the respective list entries.
    ax.set_xticklabels(col_labels, fontsize=fontsize, ha="right")
    ax.set_yticklabels(row_labels, fontsize=fontsize)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=False, bottom=True,
                   labeltop=False, labelbottom=True)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=90,
             rotation_mode="anchor")

    # Annotate
    if annotate is not None:
        if isinstance(annotate, int):
            ndecimals = annotate
        else:
            ndecimals = 3
        annotate_heatmap(array, ax, ndecimals=ndecimals, fontsize=fontsize)

    # Turn spines off and create white grid.
    for edge, spine in ax.spines.items():
        spine.set_visible(False)

    # set xticks
    ax.set_xticks(np.arange(data.shape[1] + 1) - .5, minor=True)
    ax.set_yticks(np.arange(data.shape[0] + 1) - .5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=0.5)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar, ax


def plot_variant_heatmap(arr, seq, N_res_per_heatmap_row, aaList, seq_name=None, savefig=None, figtitle=None):
    import matplotlib.pyplot as plt
    # Visualize the heatmaps
    seq_len = len(seq)
    residue_num = list(np.arange(1, seq_len + 1))
    num_heatmaps = int(np.ceil(seq_len / N_res_per_heatmap_row))
    heatmap_min = np.min(arr)
    heatmap_max = np.max(arr)
    fig, ax = plt.subplots(num_heatmaps, 1, figsize=(N_res_per_heatmap_row / len(aaList) * 4, num_heatmaps * 4))
    for k in range(num_heatmaps):
        if num_heatmaps == 1:
            ax_k = ax
        else:
            ax_k = ax[k]
        residue_num_k = residue_num[k * N_res_per_heatmap_row:min((k + 1) * N_res_per_heatmap_row, seq_len)]
        start_idx = k * N_res_per_heatmap_row
        end_idx = min((k + 1) * N_res_per_heatmap_row, seq_len)
        heatmap_k = arr[:, start_idx:end_idx]
        wt_idxs_k = np.array([[aaList.index(wt_aa),res_idx] for res_idx, wt_aa in enumerate(seq[start_idx:end_idx])])
        im = ax_k.imshow(heatmap_k, cmap="viridis", aspect="auto", vmin=heatmap_min, vmax=heatmap_max)
        ax_k.scatter(wt_idxs_k[:,1], wt_idxs_k[:,0], c='r', s=4)
        ax_k.set_yticks(range(20), aaList)
        ax_k.set_xticks(range(len(residue_num_k)), residue_num_k, fontsize=7, rotation=45)
    fig.colorbar(im, orientation='vertical')
    if figtitle is not None:
        if seq_name is not None:
            figtitle = seq_name + ': ' + figtitle
        plt.suptitle(figtitle, y=0.93, fontsize=16)
    if savefig is not None:
        plt.savefig(savefig, dpi=300, bbox_inches='tight')
    plt.show()
