# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 14:51:03 2021

@author: Jhoann
"""
import yasara
import os
from variables import address_dict, subfolders, aaList
from utils import opsys, split_mutation, mkDir, findProcess, exit_program, read_fasta

def mutate_residue(
        mutation,
        input_pdb_fname,
        log_fname,
        input_dir,
        output_dir,
        ff='AMBER15FB',
        move='all',
        mvdist=8,
        rep=1,
        surfout=0.65,
        surf=65,
        cntions=1,
        JToUnit=1.43932620865201e20, # If counterions=1, counter ions will be implicitly considered by setting the net charges to 0
        multiprocessing_proc_num=None,
):
    # format mutation / mutant
    if isinstance(mutation, list):
        mutant = '-'.join(mutation)
    else:
        mutant = mutation
        mutation = [mutation]

    # set parameters
    surfcost = (surfout) / 6.02214199e20 * JToUnit
    log_fpath = output_dir + log_fname

    # initialize lists for storing calculated values
    yasaraPDB_opt = None
    ebindDDG_opt = None

    # initialize YASARA
    print('Started YASARA')
    yasara.info.mode = 'txt'
    yasara.info.licenseshown = 0
    yasara.Console('Off')
    yasara.EnergyUnit('kcal/mol')

    if opsys == 'Windows':
        yasara_pids = findProcess('YASARA.exe')
        yasara_pid = yasara_pids[-1]
    else:
        yasara_pids = findProcess('yasara')
        yasara_pid = yasara_pids[0]
    print('All Yasara PIDs:', yasara_pids)
    print('PID for yasara program:', yasara_pid)

    yasara.Clear()

    # load input structure
    cpx = input_dir + input_pdb_fname
    yasara.LoadPDB(cpx)
    yasara.Console('Off')
    yasara.CleanAll()
    yasara.CellAuto(extension='10')
    yasara.ForceField(ff, setpar='yes')
    yasara.FixAll()

    # get mutation to WT AA and MT residue
    for mutgrp in mutation:
        WT, position, MT = split_mutation(mutgrp)
        yasara.FreeAtom(move + ' with distance< ' + str(mvdist) + ' from res ' + str(position))

    for n in range(rep):
        setname = '{0}_{1}_{2}_d{3}_s{4}_c{5}_r{6}'.format(mutant, ff, move, mvdist, surf, cntions, n)
        yasara.RandomSeed(1234 * n)

        ## WILDTYPE ##
        # swap to wildtype residue
        for mutgrp in mutation:
            WT, position, MT = split_mutation(mutgrp)
            selection = str(position)
            yasara.SwapRes(selection, WT)
            # minimize energy
            if WT != 'Gly' and WT != 'Ala': # why do we not need to OptimizeRes for Gly or Ala?
                yasara.OptimizeRes(WT + ' ' + selection, method='SCWALL')
                yasara.Boundary(Type='Periodic')
        yasara.ExperimentMinimization(convergence=0.01)
        yasara.Experiment('On')
        yasara.Wait('ExpEnd')
        yasara.ChargeObj('All', 0)
        # get energies
        epotListWT_ = yasara.Energy('All')
        epotWT_ = sum(epotListWT_)
        esolcolWT_, esolvdwWT_ = yasara.SolvEnergy(method='BoundaryFast')
        yasara.Sim('On')
        surfaccListWT_ = yasara.Surf('Accessible')
        surfaccWT_ = surfaccListWT_.pop()
        # calculate total energies
        esolWT_ = esolcolWT_ + esolvdwWT_ + surfaccWT_ * surfcost
        ebindingWT = epotWT_ + esolWT_
        yasara.Sim('Off')

        ## MUTANT ##
        # swap to mutant residue
        for mutgrp in mutation:
            WT, position, MT = split_mutation(mutgrp)
            selection = str(position)
            yasara.SwapRes(selection, MT)
            if MT != 'Gly' and MT != 'Ala':
                yasara.OptimizeRes(MT + ' ' + selection, method='SCWALL')
                yasara.Boundary(Type='Periodic')
        yasara.ExperimentMinimization(convergence=0.01)
        yasara.Experiment('On')
        yasara.Wait('ExpEnd')
        yasara.ChargeObj('All', 0)
        # get energies
        epotListMT_ = yasara.Energy('All')
        epotMT_ = sum(epotListMT_)
        yasara.Sim('On')
        esolcolMT_, esolvdwMT_ = yasara.SolvEnergy(method='BoundaryFast')
        surfaccListMT_ = yasara.Surf('Accessible')
        surfaccMT_ = surfaccListMT_.pop()
        # calculate total energies
        esolMT_ = esolcolMT_ + esolvdwMT_ + surfaccMT_ * surfcost
        ebindingMT = epotMT_ + esolMT_
        ebindingDDG = ebindingWT - ebindingMT
        # add mutated structure
        yasara.AddObj('All')
        yasara.Sim('Off')
        
        # save mutated PDB
        PDBfile_mutated = output_dir + setname + ".pdb"
        PDBname = setname + ".pdb"
        if n==0 or ebindingDDG < ebindDDG_opt:
            ebindDDG_opt = ebindingDDG
            yasaraPDB_opt = PDBname
            yasara.SavePDB(1, PDBfile_mutated)
            print("New Lowest: " + str(ebindDDG_opt))
        else:
            print("Same Lowest: " + str(ebindDDG_opt))
        yasara.LogAs(log_fpath, 'yes')
        yasara.Print('{0} {1}'.format(setname, ebindDDG_opt))

    if multiprocessing_proc_num is not None:
        print('Exiting YASARA program', yasara_pid, mutant)
        exit_program(yasara_pid)

    return yasaraPDB_opt, yasara_pid

def main():
    # define directories
    cwd = os.getcwd()
    print('Current working directory:', cwd)
    data_folder = address_dict['PIPS']
    subfolder = 'SingleMutants'
    stability_dir = data_folder+subfolders['stability']+subfolder+'/'
    pdb_dir = data_folder+subfolders['pdb']+'Input_StabilityCalc/'
    input_dir = stability_dir + 'Input/'
    outputSwapRes_dir = input_dir + 'SwapRes/'
    remove_existing_dir = False

    # define filenames
    input_fname = 'GOh1052_mutPos_DomainIII.txt'
    input_pdb_fname = 'YASARA_2EIE_GOh1052.pdb'
    log_fname = 'DDG_PKS' # 'DDG_split_2EIE'

    # get mutations
    # input is a list of positions to mutate
    if input_fname.find('.txt')>-1:
        with open(input_dir + input_fname) as f:
            wildtype_list = [mut.replace('\n','') for mut in f.readlines()]
    # input is the protein sequence --> mutate every position
    elif input_fname.find('.fasta')>-1:
        seq_dict = read_fasta(input_dir + input_fname)
        seq = list(seq_dict.values())[0]
        wildtype_list = [aa+str(i+1) for i,aa in enumerate(seq)]
    mutations = [wt+aa for wt in wildtype_list for aa in aaList if wt[0]!=aa]
    print(mutations)
    mutations = mutations[:20]
    print(mutations)

    for mutation in mutations:
        print('Running yasara_MutRes on ' + mutation + '...')
        # make SwapRes dir
        outputSwapRes_path = mkDir(mutation, outputSwapRes_dir, remove_existing_dir) + '/'
        yasaraPDB_opt, yasara_pid = mutate_residue(mutation, input_pdb_fname, log_fname, pdb_dir, outputSwapRes_path)
    print('Finished running yasara_MutRes.')



if __name__ == "__main__":  # confirms that the code is under main function
    main()


