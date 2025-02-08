# -*- coding: utf-8 -*-
import yasara
import numpy as np
import os
import time
from variables import address_dict, subfolders, aaList
from utils import opsys, split_mutation, findProcess, exit_program, save_dict_as_csv, combine_csv_files


# CALCULATE BINDING ENERGY INCLUDING IMPLICIT SOLVATION
# Created by Sebastian Maurer-Stroh, inspired by template macros from Elmar Krieger.
# Edited by Jhoann Miyajima

def get_mutstr(mutation):
    if isinstance(mutation, list):
        mutstr = '+'.join(mutation)
    else:
        mutstr = mutation
        mutation = [mutation]
    return mutstr, mutation
def get_res_to_mutate_from_inputfile(col, input_fname, input_dir):
    res_to_mutate = []
    if input_fname is not None:
        import pandas as pd
        # get positions to mutate from input file
        df = pd.read_csv(input_dir + input_fname)
        res_to_mutate = df[col].dropna().tolist()
    return res_to_mutate
def get_mutations_list(res_to_mutate, aaList=aaList):
    # check if amino acid to mutate to is already specified
    if res_to_mutate[0][-1] in aaList:
        mutations_list = res_to_mutate
    else:
        mutations_list = []
        for wt in res_to_mutate:
            wtAA = wt[0]  # Get the first character of the sequence
            for aa in aaList:
                if aa != wtAA:
                    mt = wt + aa
                    mutations_list.append(mt)
    print(len(mutations_list), mutations_list)
    return mutations_list

def set_up_sce_for_minimization(mutation, move, mvdist, mvdrug, continue_processing=True):
    # allow mutated positions to move
    mutname = ''
    for mutgrp in mutation:
        WT, position, MT = split_mutation(mutgrp)
        # check if current residue position is the same
        WT_actual = yasara.NameRes(position)[0]
        print('Actual WT residue:', WT_actual, '; Target WT residue:', WT.upper())
        if WT_actual != WT.upper():
            print('Incorrect WT amino acid found at position to mutate >> End processing')
            continue_processing = False
            break
        else:
            print('Correct WT amino acid found at position to mutate >> Continue processing')
        x = '{0}{1}{2}'.format(WT, position, MT)
        mutname += x
        yasara.ShowRes(str(position))
        yasara.FreeAtom(move + ' with distance <' + str(mvdist) +' from res ' + str(position))
    print('mutname:', mutname)

    if (mvdrug == 1): # allow ligand to move
            yasara.FreeAtom(move + ' with distance <' + str(mvdist) +' from Obj 2')
    if (mvdrug == 0): # fix ligand
            yasara.FixObj("2")
    return mutname, continue_processing

def save_sce_as_pdb(output_dir, struct_subdir='postOpt/', del_sce_files=False):
    # convert sce files to pdb files
    sce_list = [f for f in os.listdir(output_dir + struct_subdir) if f.endswith('.sce')]
    for sce_fname in sce_list:
        yasara.Clear()
        yasara.LoadSce(output_dir + struct_subdir + sce_fname)
        yasara.JoinObj('2', '1', center='No')
        yasara.NameMol('Protein', 'R')
        yasara.NameMol('not Protein', 'L')
        yasara.SavePDB('1', output_dir + struct_subdir + sce_fname.replace('.sce', '.pdb'))
        print('Converted sce to pdb file:', sce_fname)
    if del_sce_files:
        for sce_fname in sce_list:
            os.remove(output_dir+struct_subdir+sce_fname)

def mutate_residue(
    mutation,
    struct_fname,
    struct_dir,
    output_dir,
    move='!backbone',
    minimize_energy=True,
    resetSce=False,
    multiprocessing_proc_num=None,
    nrep=5,
    ff = 'AMBER15FB',
    mvdist = 4,
    mvdrug = 1,
    surfout = 1.75,
    cntions = 0,
    JToUnit = 1.43932620865201e20,  # If counterions=1, counter ions will be implicitly considered by setting the net charges to 0
    structure_dict = {'Cpx': 3, 'Rtr': 1, 'Lgd': 2},
    yasara_pid=None,
):

    # format mutation / mutant
    mutant, mutation = get_mutstr(mutation)

    # settings and fpaths
    surfcost = (surfout)/6.02214199e20*JToUnit
    lig = struct_fname.split("_")[0]
    ligname = lig
    struct = struct_fname.split("_")[2]
    struct_fpath = struct_dir + struct_fname + '.sce'
    log_fname = 'DDG_' + struct
    log_fpath = output_dir + log_fname
    element_list = ['Rtr', 'Lgd', 'Cpx']
    setname = '{0}_{1}_{2}_{3}_rs{4}_{5}_d{6}_md{7}_s{8}_c{9}_r{10}'.format(ligname, struct, mutant, ff,
                                                                                   resetSce * 1, move, mvdist, mvdrug,
                                                                                   surfout, cntions, nrep)

    # INITIALIZE ENERGY FEATURES RESULTS DICT #
    energy_features_list = ['ebindDDG', 'ebindMT', 'ebindWT'] + \
                            [feature+obj for feature in ['ebindWT', 'epotWT', 'esolWT', 'esolWTcol', 'esolWTvdw', 'surfaccWT'] for obj in element_list] + \
                            [feature+obj for feature in ['ebindMT', 'epotMT', 'esolMT', 'esolMTcol', 'esolMTvdw', 'surfaccMT'] for obj in element_list]
    key_energy_features = ['ebindDDG', 'ebindMT', 'ebindWT', 'ebindMTCpx', 'ebindMTRtr', 'ebindMTLgd', 'ebindWTCpx', 'ebindWTRtr', 'ebindWTLgd']
    res_avg_cols = ['struct', 'mutant', 'ligname', 'setname', 'minimize_energy', 'resetSce', 'move', 'mvdist'] + \
                key_energy_features + [f+'_std' for f in key_energy_features] + \
                [feature + obj for feature in ['epotWT', 'esolWT', 'esolWTcol', 'esolWTvdw', 'surfaccWT'] for obj in element_list] + \
                [feature + obj for feature in ['epotMT', 'esolMT', 'esolMTcol', 'esolMTvdw', 'surfaccMT'] for obj in element_list] + \
                ['surfcost', 'mvdrug', 'ff', 'counterions', 'nrep']
    res_all_cols = ['struct', 'mutant', 'ligname', 'setname', 'n', 'surfcost', 'minimize_energy', 'resetSce', 'move', 'mvdist', 'mvdrug', 'ff', 'counterions'] + energy_features_list

    res = {}
    for f in energy_features_list: res.update({f:[]})
    res_avg = {'struct':struct, 'mutant':mutant, 'ligname':ligname, 'setname':setname,
               'surfcost':surfcost, 'minimize_energy': minimize_energy, 'resetSce':resetSce, 'move':move, 'mvdist':mvdist, 'mvdrug':mvdrug, 'ff':ff, 'counterions':cntions, 'nrep':nrep}

    # INITIALIZE YASARA #
    print('Starting Yasara processing for ' + struct_fname + ' >> ' + mutant)
    yasara.info.mode = 'txt'
    print('Set info mode')
    yasara.info.licenseshown = 0
    print('Set license shown')
    time = yasara.SystemTime()
    print('Set system time')
    # get yasara PID
    input_yasara_pid = yasara_pid
    if input_yasara_pid is None:
        if opsys=='Windows':
            yasara_pids = findProcess('YASARA.exe')
            yasara_pid = yasara_pids[-1]
        else:
            yasara_pids = findProcess('yasara')
            yasara_pid = yasara_pids[0]
        print('All Yasara PIDs:', yasara_pids)
        print('PID for yasara program:', yasara_pid)
    yasara.info.pid = int(time)+5
    print('Set info pid')
    yasara.Console('Off')
    yasara.Processors(1)
    print('Set Processors')
    yasara.EnergyUnit('kcal/mol')
    print('Set Energy Unit')
    yasara.ForceField(ff, setpar='yes')
    print('Set Force Field')

    # allow mutated positions to move
    continue_processing = True
    # iterate through energy minimization runs
    for n in range(nrep):
        print('rep#='+str(n))
        setname_n = setname + '-' + str(n)
        # clear and load scene
        yasara.RandomSeed(1234*n)
        if n==0 or resetSce:
            yasara.Clear()
            print('Clear previous contents')
            yasara.LoadSce(struct_fpath)
            print('Load Sce')
            yasara.CleanAll()
            print('Clean All')
            yasara.CellAuto(extension='10')
            yasara.FixAll()
            # set up scene for minimization
            mutname, continue_processing = set_up_sce_for_minimization(mutation, move, mvdist, mvdrug, continue_processing)

        if not continue_processing:
            print('File could not be processed. Moving on to the next simulation...')
            break
        else:
            # perform energy minimization calculations for WT and MT
            for WT_or_MT in ['WT', 'MT']:

                # perform mutagenesis
                print('Performing SwapRes for ' + WT_or_MT + '...')
                for mutgrp in mutation:
                    print(mutgrp)
                    WT, position, MT = split_mutation(mutgrp)
                    selection = str(position)
                    if WT_or_MT=='MT':
                        aa_for_swapres = MT
                    elif WT_or_MT=='WT':
                        aa_for_swapres = WT
                    yasara.SwapRes(selection, aa_for_swapres)
                    if aa_for_swapres!='Gly' and aa_for_swapres!='Ala':
                        yasara.OptimizeRes(aa_for_swapres + ' ' + selection, method = 'SCWALL')
                        yasara.Boundary(Type='Periodic')

                # perform energy minimization
                if minimize_energy:
                    print('Performing Energy Minimization for ' + WT_or_MT + '...')
                    yasara.ExperimentMinimization(convergence=0.01)
                    yasara.Experiment('On')
                    yasara.Wait('ExpEnd')
                    print('ExpEnd')

                print('Calculating energies for ' + WT_or_MT + '...')
                for element, element_idx in structure_dict.items():
                    print('Processing element ' + str(element_idx) + ': ' + element)
                    # remove all except that element
                    if element in ['Rtr', 'Lgd']:
                        yasara.RemoveObj('not ' + str(element_idx))
                    # set charge to 0
                    yasara.ChargeObj('All',0)
                    # get potential energy
                    epotList_ = yasara.Energy('All')
                    epot_ = sum(epotList_)
                    # get component solvation energies
                    esolcol_,esolvdw_ = yasara.SolvEnergy(method='BoundaryFast')
                    yasara.Sim('On')
                    # get interfacial energy
                    surfaccList_ = yasara.Surf('Accessible')
                    surfacc_ = surfaccList_.pop()
                    # get total solvation energy
                    esol_ = esolcol_+esolvdw_+surfacc_*surfcost
                    yasara.AddObj('All')

                    # update results
                    res['ebind' + WT_or_MT + element].append(epot_ + esol_)
                    res['epot' + WT_or_MT + element].append(epot_)
                    res['esol' + WT_or_MT + element].append(esol_)
                    res['esol' + WT_or_MT + 'col' + element].append(esolcol_)
                    res['esol' + WT_or_MT + 'vdw' + element].append(esolvdw_)
                    res['surfacc' + WT_or_MT + element].append(surfacc_)

                # DG_bind = DG_complex - (DG_receptor + DG_ligand)
                res['ebind' + WT_or_MT].append(res['ebind'+WT_or_MT+'Cpx'][-1] - res['ebind'+WT_or_MT+'Rtr'][-1] - res['ebind'+WT_or_MT+'Lgd'][-1])
                yasara.Sim('Off')
                print('Obtained energies for ' + WT_or_MT)
                print('ebinding' + WT_or_MT + ':', res['ebind'+WT_or_MT])

                # save structure as sce
                setname_mod = setname_n[:setname_n.find('_')] + '_postOpt-'+WT_or_MT + setname_n[setname_n.find('_'):]
                output_postOpt_fpath = output_dir + 'postOpt/' + setname_mod
                yasara.SaveSce(output_postOpt_fpath + '.sce')

            # calculate energy change
            res['ebindDDG'].append(res['ebindMT'][-1] - res['ebindWT'][-1])
            yasara.Sim('Off')


    if continue_processing:
        # DDG
        # Get average and standard deviation for energies calculated
        for f in energy_features_list:
            res_avg[f] = round(np.mean(np.array(res[f])),4)
            if f in key_energy_features:
                res_avg[f+'_std'] = round(np.std(np.array(res[f])),4)
        print("Obtained average energy changes")
        print('ebindDDG=' + str(round(res_avg['ebindDDG'],2)))

        # save average results as csv
        csv_txt_avg, log_fpath_avg, write_mode = save_dict_as_csv(res_avg, res_avg_cols, log_fpath, multiprocessing_proc_num=multiprocessing_proc_num)
        print('Saved AVG results for ' + struct + mutant + ' to CSV (mode=' + write_mode + ').')

        # save all results as csv
        res.update({'struct': [struct]*nrep, 'mutant': [mutant]*nrep, 'ligname': [ligname]*nrep, 'setname': [setname+'-'+str(n) for n in range(nrep)], 'n': [n for n in range(nrep)],
                     'surfcost': [surfcost]*nrep, 'minimize_energy': [minimize_energy]*nrep, 'resetSce': [resetSce]*nrep, 'move': [move]*nrep,
                     'mvdist': [mvdist]*nrep, 'mvdrug': [mvdrug]*nrep, 'ff': [ff]*nrep, 'counterions': [cntions]*nrep})
        csv_txt_avg, log_fpath_full, write_mode = save_dict_as_csv(res, res_all_cols, log_fpath, csv_suffix='_FULL', multiprocessing_proc_num=multiprocessing_proc_num)
        print('Saved FULL results for ' + struct + mutant + ' to CSV (mode=' + write_mode + ').')

    if multiprocessing_proc_num is not None:
        # convert sce to pdb files
        save_sce_as_pdb(output_dir, struct_subdir='postOpt/', del_sce_files=True)
        # exit program
        print('Exiting YASARA program')
        exit_program(yasara_pid)

    return yasara_pid, log_fpath_avg, log_fpath_full


def main():

    # set inputs and parameters
    run_multiprocessing =  15 # None #
    nrep = 10 # 2 #
    data_folder = address_dict['influenza-resistance'] # address_dict['PIPS']
    subfolder = 'PA-A37T' # 'NA' #
    output_dir = data_folder + subfolders['yasara'] + 'Output/' + subfolder + '/' # data_folder + subfolders['yasara'] + 'Output/DomainII_output/'
    input_dir = data_folder + subfolders['yasara'] + 'Input/'
    struct_dir = data_folder + subfolders["sce"] + subfolder + '/'
    output_fname = 'DDG_PA-A37T_Baloxavir' # 'DDG_NA_Oseltamivir' #  None
    input_fname = None # "HTP_DomainPos_Final.csv"
    input_col = 'Domain III'
    if input_fname is not None:
        res_to_mutate = get_res_to_mutate_from_inputfile(input_col, input_fname, input_dir)
    struct_to_mutate_dict = {
        # NA
        # 'Oseltamivir_preOpt_NA-H1N1-3TI6-1': ['V263T', 'H274Y', ['V263T', 'H274Y']],
        # 'Oseltamivir_preOpt_NA-H1N1-3TI6-2': ['V263T', 'H274Y', ['V263T', 'H274Y']],
        # 'Oseltamivir_preOpt_NA-H1N1-Victoria1162A-5NWE-1': ['T264V', 'Y275H', ['T264V', 'Y275H']],
        # 'Oseltamivir_preOpt_NA-H1N1-Victoria1162A-5NWE-2': ['T264V', 'Y275H', ['T264V', 'Y275H']],
        # 'Oseltamivir_preOpt_NA-H1N1-Victoria1162A-5NWE-3': ['T264V', 'Y275H', ['T264V', 'Y275H']],
        # 'Oseltamivir_preOpt_NA-H1N1-Victoria1162A-5NWE-4': ['T264V', 'Y275H', ['T264V', 'Y275H']],
        # 'Oseltamivir_preOpt_NA-H1N1-Wisconsin67-5NWE-1': ['I264T', 'I264V', 'H275Y', ['I264T', 'H275Y'], ['I264V', 'H275Y']],
        # 'Oseltamivir_preOpt_NA-H1N1-Wisconsin67-5NWE-2': ['I264T', 'I264V', 'H275Y', ['I264T', 'H275Y'], ['I264V', 'H275Y']],
        # 'Oseltamivir_preOpt_NA-H1N1-Wisconsin67-5NWE-3': ['I264T', 'I264V', 'H275Y', ['I264T', 'H275Y'], ['I264V', 'H275Y']],
        # 'Oseltamivir_preOpt_NA-H1N1-Wisconsin67-5NWE-4': ['I264T', 'I264V', 'H275Y', ['I264T', 'H275Y'], ['I264V', 'H275Y']],
        # 'Oseltamivir_preOpt_NA-H1N1-Victoria1162A-5NZN': ['V264T', 'Y275H', ['V264T', 'Y275H']],
        # 'Oseltamivir_preOpt_NA-H1N1-Wisconsin67-5NZN': ['V264T', 'Y275H', ['V264T', 'Y275H']],

        # PA
        'Baloxavir_preOpt_PA-H3N2-res1_fixed': ['A37T'],
        'Baloxavir_preOpt_PA-H3N2-Croatia10136RV-8T5V': ['A37T'],
        'Baloxavir_preOpt_PA-H3N2-Tasmania104-8T5V-Croatia10136RVmodel': ['T37A']
    }

    # run conditions
    minimize_energy_list = [True]
    resetSce_list = [False, True]
    move_list = ['!backbone', 'all']

    # create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print('Created output directory:', output_dir)
    if not os.path.exists(output_dir+'postOpt/'):
        os.makedirs(output_dir+'postOpt/')
        print('Created postOpt sub-directory:', output_dir+'postOpt/')

    if run_multiprocessing is None:
        proc_num = 0
        current_yasara_pid = None
        log_fpath_list = []
        log_fpath_list_FULL = []
        for struct_fname, res_to_mutate in struct_to_mutate_dict.items():
            mutations_list = get_mutations_list(res_to_mutate, aaList=aaList)
            # process structure with mutations list
            for resetSce in resetSce_list:
                for minimize_energy in minimize_energy_list:
                    for move in move_list:
                        for i, mutation in enumerate(mutations_list):
                            mutstr, mutation = get_mutstr(mutation)
                            print('Processing Struct ' + struct_fname + 'Mutation # ' + str(i) + ': ' + mutstr + '; Move: ' + move + '; Minimize energy: ' + str(minimize_energy))
                            current_yasara_pid, log_fpath_avg, log_fpath_full = mutate_residue(mutation, struct_fname, struct_dir, output_dir, move,
                                                                minimize_energy, resetSce, yasara_pid=current_yasara_pid, nrep=nrep)
                            if log_fpath_avg not in log_fpath_list:
                                log_fpath_list.append(log_fpath_avg)
                            if log_fpath_full not in log_fpath_list_FULL:
                                log_fpath_list_FULL.append(log_fpath_full)
                            proc_num += 1

        # convert sce to pdb files
        save_sce_as_pdb(output_dir, struct_subdir='postOpt/', del_sce_files=True)
        # exit program
        exit_program(current_yasara_pid)

    else:
        from multiprocessing import Pool, Process, cpu_count
        # get list of conditions
        args_list = []
        log_fpath_list = []
        log_fpath_list_FULL = []
        proc_num = 0
        for struct_fname, res_to_mutate in struct_to_mutate_dict.items():
            mutations_list = get_mutations_list(res_to_mutate, aaList=aaList)
            for resetSce in resetSce_list:
                for minimize_energy in minimize_energy_list:
                    for move in move_list:
                        for mutation in mutations_list:
                            args_list.append((mutation, struct_fname, struct_dir, output_dir, move, minimize_energy, resetSce, proc_num, nrep))
                            log_fpath_avg = output_dir + 'DDG_' + struct_fname.split("_")[2] + '_' + str(proc_num) + '.csv'
                            log_fpath_full = output_dir + 'DDG_' + struct_fname.split("_")[2] + '_FULL_' + str(proc_num) + '.csv'
                            if log_fpath_avg not in log_fpath_list:
                                log_fpath_list.append(log_fpath_avg)
                            if log_fpath_full not in log_fpath_list_FULL:
                                log_fpath_list_FULL.append(log_fpath_full)
                            proc_num += 1

        # MULTIPROCESSING POOL APPLY_ASYNC
        num_cpu = cpu_count()
        print("Number of cpu : ", num_cpu)
        pool = Pool(processes=num_cpu)
        for proc_num, args in enumerate(args_list):
            (mutation, struct_fname, struct_dir, output_dir, move, minimize_energy, resetSce, proc_num, nrep) = args
            mutstr, mutation = get_mutstr(mutation)
            pool.apply_async(mutate_residue, args=args)
            time.sleep(3)
            print('Spawned process # ' + str(proc_num) + '. Mutation:' + mutstr + '; Move: ' + move + '; Minimize energy: ' + str(minimize_energy))
        pool.close()
        try:
            pool.join()
        except:
            print('Unable to join pool. Moving on...')
        print('Completed all ' + str(proc_num) + 'processes.')

    # combine files spawned
    ## AVG results
    _ = combine_csv_files(log_fpath_list, output_dir, output_fname)
    print('Saved all AVG results as CSV:' + output_dir + output_fname + '.csv')
    
    ## FULL results
    _ = combine_csv_files(log_fpath_list_FULL, output_dir, output_fname+'_FULL')
    print('Saved all FULL results as CSV:' + output_dir + output_fname + '_FULL' + '.csv')

if __name__ == "__main__":  # confirms that the code is under main function
    main()
