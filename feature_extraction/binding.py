import yasara
import numpy as np
import os
if os.path.basename(os.getcwd()) != 'feature_extraction': os.chdir('./feature_extraction/')
from utils.utils import aaList, mkDir, opsys, split_mutation, get_mutstr, findProcess, exit_program, save_dict_as_csv, combine_csv_files
if opsys == 'Windows': yasara_process_name = 'YASARA.exe'
else: yasara_process_name = 'yasara'


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
    # else specify all the possible mutations to evaluate
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
        yasara.FreeAtom(move + ' with distance <' + str(mvdist) + ' from res ' + str(position))
    print('mutname:', mutname)

    if (mvdrug == 1):  # allow ligand to move
        yasara.FreeAtom(move + ' with distance <' + str(mvdist) + ' from Obj 2')
    if (mvdrug == 0):  # fix ligand
        yasara.FixObj("2")
    return mutname, continue_processing

def mutate_residue(
        mutation,
        struct_fname,
        struct_dir,
        output_dir,
        move='!backbone',
        minimize_energy=True,
        resetSce=False,
        nrep=5,
        ff='AMBER15FB',
        mvdist=4,
        mvdrug=1,
        surfout=1.75,
        cntions=0,
        JToUnit=1.43932620865201e20,
        # If counterions=1, counter ions will be implicitly considered by setting the net charges to 0
        structure_dict={'Cpx': 3, 'Rtr': 1, 'Lgd': 2},
):
    # format mutation / mutant
    mutant, mutation = get_mutstr(mutation)

    # settings and fpaths
    surfcost = (surfout) / 6.02214199e20 * JToUnit
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
                           [feature + obj for feature in
                            ['ebindWT', 'epotWT', 'esolWT', 'esolWTcol', 'esolWTvdw', 'surfaccWT'] for obj in
                            element_list] + \
                           [feature + obj for feature in
                            ['ebindMT', 'epotMT', 'esolMT', 'esolMTcol', 'esolMTvdw', 'surfaccMT'] for obj in
                            element_list]
    key_energy_features = ['ebindDDG', 'ebindMT', 'ebindWT', 'ebindMTCpx', 'ebindMTRtr', 'ebindMTLgd', 'ebindWTCpx',
                           'ebindWTRtr', 'ebindWTLgd']
    res_avg_cols = ['struct', 'mutant', 'ligname', 'setname', 'minimize_energy', 'resetSce', 'move', 'mvdist'] + \
                   key_energy_features + [f + '_std' for f in key_energy_features] + \
                   [feature + obj for feature in ['epotWT', 'esolWT', 'esolWTcol', 'esolWTvdw', 'surfaccWT'] for obj in
                    element_list] + \
                   [feature + obj for feature in ['epotMT', 'esolMT', 'esolMTcol', 'esolMTvdw', 'surfaccMT'] for obj in
                    element_list] + \
                   ['surfcost', 'mvdrug', 'ff', 'counterions', 'nrep']
    res_all_cols = ['struct', 'mutant', 'ligname', 'setname', 'n', 'surfcost', 'minimize_energy', 'resetSce', 'move',
                    'mvdist', 'mvdrug', 'ff', 'counterions'] + energy_features_list

    res = {}
    for f in energy_features_list: res.update({f: []})
    res_avg = {'struct': struct, 'mutant': mutant, 'ligname': ligname, 'setname': setname,
               'surfcost': surfcost, 'minimize_energy': minimize_energy, 'resetSce': resetSce, 'move': move,
               'mvdist': mvdist, 'mvdrug': mvdrug, 'ff': ff, 'counterions': cntions, 'nrep': nrep}

    # INITIALIZE YASARA #
    print('Starting Yasara processing for ' + struct_fname + ' >> ' + mutant)
    yasara.info.mode = 'txt'
    print('Set info mode')
    yasara.info.licenseshown = 0
    print('Set license shown')
    time = yasara.SystemTime()
    print('Set system time')
    yasara.info.pid = int(time) + 5
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
        print('rep#=' + str(n))
        setname_n = setname + '-' + str(n)
        # clear and load scene
        yasara.RandomSeed(1234 * n)
        if n == 0 or resetSce:
            yasara.Clear()
            print('Clear previous contents')
            yasara.LoadSce(struct_fpath)
            print('Load Sce')
            yasara.CleanAll()
            print('Clean All')
            yasara.CellAuto(extension='10')
            yasara.FixAll()
            # set up scene for minimization
            mutname, continue_processing = set_up_sce_for_minimization(mutation, move, mvdist, mvdrug,
                                                                       continue_processing)

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
                    if WT_or_MT == 'MT':
                        aa_for_swapres = MT
                    elif WT_or_MT == 'WT':
                        aa_for_swapres = WT
                    yasara.SwapRes(selection, aa_for_swapres)
                    if aa_for_swapres != 'Gly' and aa_for_swapres != 'Ala':
                        yasara.OptimizeRes(aa_for_swapres + ' ' + selection, method='SCWALL')
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
                    yasara.ChargeObj('All', 0)
                    # get potential energy
                    epotList_ = yasara.Energy('All')
                    epot_ = sum(epotList_)
                    # get component solvation energies
                    esolcol_, esolvdw_ = yasara.SolvEnergy(method='BoundaryFast')
                    yasara.Sim('On')
                    # get interfacial energy
                    surfaccList_ = yasara.Surf('Accessible')
                    surfacc_ = surfaccList_.pop()
                    # get total solvation energy
                    esol_ = esolcol_ + esolvdw_ + surfacc_ * surfcost
                    yasara.AddObj('All')

                    # update results
                    res['ebind' + WT_or_MT + element].append(epot_ + esol_)
                    res['epot' + WT_or_MT + element].append(epot_)
                    res['esol' + WT_or_MT + element].append(esol_)
                    res['esol' + WT_or_MT + 'col' + element].append(esolcol_)
                    res['esol' + WT_or_MT + 'vdw' + element].append(esolvdw_)
                    res['surfacc' + WT_or_MT + element].append(surfacc_)

                # DG_bind = DG_complex - (DG_receptor + DG_ligand)
                res['ebind' + WT_or_MT].append(
                    res['ebind' + WT_or_MT + 'Cpx'][-1] - res['ebind' + WT_or_MT + 'Rtr'][-1] -
                    res['ebind' + WT_or_MT + 'Lgd'][-1])
                yasara.Sim('Off')
                print('Obtained energies for ' + WT_or_MT)
                print('ebinding' + WT_or_MT + ':', res['ebind' + WT_or_MT])

                # save structure as sce
                setname_mod = setname_n[:setname_n.find('_')] + '_postOpt-' + WT_or_MT + setname_n[setname_n.find('_'):]
                output_postOpt_fpath = output_dir + 'postOpt/' + setname_mod
                yasara.SaveSce(output_postOpt_fpath + '.sce')

            # calculate energy change
            res['ebindDDG'].append(res['ebindMT'][-1] - res['ebindWT'][-1])
            yasara.Sim('Off')

    if continue_processing:
        # DDG
        # Get average and standard deviation for energies calculated
        for f in energy_features_list:
            res_avg[f] = round(np.mean(np.array(res[f])), 4)
            if f in key_energy_features:
                res_avg[f + '_std'] = round(np.std(np.array(res[f])), 4)
        print("Obtained average energy changes")
        print('ebindDDG=' + str(round(res_avg['ebindDDG'], 2)))

        # save average results as csv
        csv_txt_avg, log_fpath_avg, write_mode = save_dict_as_csv(res_avg, res_avg_cols, log_fpath)
        print('Saved AVG results for ' + struct + mutant + ' to CSV (mode=' + write_mode + ').')

        # save all results as csv
        res.update({'struct': [struct] * nrep, 'mutant': [mutant] * nrep, 'ligname': [ligname] * nrep,
                    'setname': [setname + '-' + str(n) for n in range(nrep)], 'n': [n for n in range(nrep)],
                    'surfcost': [surfcost] * nrep, 'minimize_energy': [minimize_energy] * nrep,
                    'resetSce': [resetSce] * nrep, 'move': [move] * nrep,
                    'mvdist': [mvdist] * nrep, 'mvdrug': [mvdrug] * nrep, 'ff': [ff] * nrep,
                    'counterions': [cntions] * nrep})
        csv_txt_avg, log_fpath_full, write_mode = save_dict_as_csv(res, res_all_cols, log_fpath, csv_suffix='_FULL')
        print('Saved FULL results for ' + struct + mutant + ' to CSV (mode=' + write_mode + ').')

    return log_fpath_avg, log_fpath_full


def get_yasara_binding_features():
    # set inputs and parameters
    nrep = 5
    input_dir = '../data/feature_extraction/Input/'
    struct_dir = input_dir
    output_dir = '../data/feature_extraction/'
    input_fname = 'HTP_DomainPos_Final.csv'
    input_col = 'Domain III'
    output_fname = 'DDGbinding_GOh1052-DomainIIImut_S82'
    struct_to_mutate_dict = {
        'S152_1GOG_GOh1001b_postOpt': get_res_to_mutate_from_inputfile(input_col, input_fname, input_dir)
    }
    mkDir('postOpt', output_dir)
    proc_num = 0
    log_fpath_list = []
    log_fpath_list_FULL = []
    for struct_fname, res_to_mutate in struct_to_mutate_dict.items():
        mutations_list = get_mutations_list(res_to_mutate, aaList=aaList)
        for i, mutation in enumerate(mutations_list):
            mutstr, mutation = get_mutstr(mutation)
            log_fpath_avg, log_fpath_full = mutate_residue(mutation, struct_fname, struct_dir, output_dir,
                                                           move='!backbone', minimize_energy=True, resetSce=False, nrep=nrep)
            if log_fpath_avg not in log_fpath_list:
                log_fpath_list.append(log_fpath_avg)
            if log_fpath_full not in log_fpath_list_FULL:
                log_fpath_list_FULL.append(log_fpath_full)
            proc_num += 1
    yasara_pid = findProcess(yasara_process_name)[-1]
    exit_program(yasara_pid)

    # combine files spawned
    _ = combine_csv_files(log_fpath_list, output_dir, output_fname)
    print('Saved all AVG results as CSV:' + output_dir + output_fname + '.csv')
    _ = combine_csv_files(log_fpath_list_FULL, output_dir, output_fname + '_FULL')
    print('Saved all FULL results as CSV:' + output_dir + output_fname + '_FULL' + '.csv')

if __name__ == "__main__":  # confirms that the code is under main function
    get_yasara_binding_features()
