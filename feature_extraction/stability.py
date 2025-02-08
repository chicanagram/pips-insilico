import os
import time
import numpy as np
from multiprocessing import Pool, cpu_count
from variables import address_dict, subfolders, aaList
from utils import opsys, mkDir, sort_list, get_mutstr, read_fasta, exit_program


def run_stability_pipeline(mutation, input_pdb_fname, process_steps, pdb_dir, swapRes_dir, repair_dir, build_dir, log_fname, receptor_molname='R', remove_existing_dir=True, foldx_abspath=os.path.abspath('../../../yasara/foldx_2025/foldx_20251231_linux'), multiprocessing_proc_num=None):
    # RUN YASARA MUTAGENESIS ##
    yasara_pid = None
    mutstr, mutation = get_mutstr(mutation)

    #  make directories
    SwapRes_path = mkDir(mutstr, swapRes_dir, remove_existing_dir) + '/'
    Repair_path = mkDir(mutstr, repair_dir, remove_existing_dir) + '/'
    Build_path = mkDir(mutstr, build_dir, remove_existing_dir) + '/'
    print('Save folders:', SwapRes_path, Repair_path, Build_path)

    if 'yasara_MutRes' in process_steps:
        from run_yasara_mutagenesis_stability import mutate_residue
        print('Running yasara_MutRes on ' + mutstr + '...')
        start_time = time.time()
        # perform YASARA mutation of residue
        print('YASARA MutRes Input PDB:', input_pdb_fname)
        swapRes_output_pdb, yasara_pid = mutate_residue(mutation, input_pdb_fname, log_fname, pdb_dir, SwapRes_path, multiprocessing_proc_num=multiprocessing_proc_num)
        end_time = time.time()
        print('Finished running YASARA mutagenesis on '+mutstr, '>> Time taken:', str(round((end_time-start_time)/60,2)) + ' min')
        print('PDB output:', swapRes_output_pdb)

    ## RUN FOLDX REPAIR ##
    if 'foldx_Repair' in process_steps:
        from stability.foldx_Repair_MutateResidue import foldxRepair
        print('Running foldxRepair on ' + ','.join(mutstr) + '...')
        start_time = time.time()
        # get PDB input
        foldx_repair_inputPDB = sort_list([f for f in os.listdir(SwapRes_path) if f.endswith('.pdb')])[0]
        print('FoldX Repair Input PDB:', foldx_repair_inputPDB)
        # perform FoldX PDB Repair
        foldxRepair(foldx_repair_inputPDB, os.path.abspath(Repair_path), os.path.abspath(SwapRes_path), foldx_abspath)
        end_time = time.time()
        print('Finished running foldxRepair on '+mutstr, '>> Time taken:', str(round((end_time-start_time)/60,2)) + ' min')

    ## RUN FOLDX BUILD ##
    if 'foldx_Build' in process_steps:
        from stability.foldx_Repair_MutateResidue import foldxBuild
        print('Running foldxBuild on ' + mutstr + '...')
        start_time = time.time()
        # get PDB input
        foldx_build_inputPDB = sort_list([f for f in os.listdir(Repair_path) if f.endswith('.pdb')])[0]
        print('FoldX Build Input PDB:', foldx_build_inputPDB)
        # make PDB folder in Build dir
        # pdbName = foldx_build_inputPDB.replace('_Repair.pdb', '')
        # pdbBuildPath = mkDir(pdbName, Build_path) + '/'
        # perform FoldX PDB Build
        # foldxBuild(mutation, pdbName, os.path.abspath(pdbBuildPath), os.path.abspath(Repair_path), foldx_abspath, receptor_molname)
        foldxBuild(mutation, foldx_build_inputPDB.strip('_Repair.pdb'), os.path.abspath(Build_path), os.path.abspath(Repair_path), foldx_abspath, receptor_molname)
        end_time = time.time()
        print('Finished running foldxBuild on '+mutstr, '>> Time taken:', str(round((end_time-start_time)/60,2)) + ' min')

    return yasara_pid

def stability_pipeline_wrapper(mutations_list, input_pdb_fname, process_steps, pdb_dir, swapRes_dir, repair_dir, build_dir, log_fname, receptor_molname, remove_existing_dir=True, foldx_abspath=os.path.abspath('../../../yasara/foldx_2025/foldx_20251231_linux'), multiprocessing_proc_num=None):
    for i, mutation in enumerate(mutations_list):
        yasara_pid_i = run_stability_pipeline(mutation, input_pdb_fname, process_steps,
                                            pdb_dir, swapRes_dir, repair_dir, build_dir, log_fname,
                                            remove_existing_dir=remove_existing_dir, foldx_abspath=foldx_abspath,
                                            multiprocessing_proc_num=None)
        print('**** Processed ' + mutation, str(i + 1) + '/' + str(len(mutations_list)) + ' ****')
        if i==0: yasara_pid = yasara_pid_i
    if yasara_pid is not None: exit_program(yasara_pid)
    return yasara_pid

def main():
    # define directories & filenames
    cwd = os.getcwd()
    print('Current working directory:', cwd)
    data_folder = address_dict['influenza-resistance'] # address_dict['SoluProtMut'] # address_dict['PIPS2']
    subfolder = 'NA_enzymeonly' # 'PKS' # 'ET096'
    input_fname = 'NA_evalmut_V264T.txt' # 'soluprotmutdb_2level_Type III PKS pyrrolidine ketide synthase.txt' # 'ETS83096.fasta' # 'GOh1052_mutPos_DomainIII.txt'
    input_pdb_fname = 'NA-H1N1-Victoria1162A-5NWE-2-264V-275Y.pdb' # 'PKS_UniProtA0A3G4RHW3_AlphaFold.pdb' # 'ETS83096.pdb' # 'YASARA_2EIE_GOh1052.pdb'
    output_fname = input_pdb_fname.strip('.pdb') + '_stability_features.csv'
    log_fname = 'DDG_' + subfolder # 'DDG_split_2EIE'
    stability_dir = data_folder + subfolders['stability'] + subfolder + '/'
    pdb_dir = data_folder + subfolders['pdb'] + subfolder + '/' # 'UPO_batch1_enzymeOnly/' # data_folder + subfolders['pdb'] + '2EIE_ProteinOnly/'
    input_dir = stability_dir + 'Input/'
    output_dir = stability_dir + 'Output/' + input_pdb_fname.strip('.pdb') + '/'
    swapRes_dir = output_dir + 'SwapRes/'
    repair_dir = output_dir + 'Repair/'
    build_dir = output_dir + 'Build/'
    remove_existing_dir = False
    run_multiprocessing = None # 15
    num_res_per_group = 2 # len(wildtype_list) # 40
    receptor_molname = 'R'

    # define process steps
    process_steps = ['yasara_MutRes', 'foldx_Repair', 'foldx_Build', 'combine_results']
    # process_steps = ['yasara_MutRes']
    # process_steps = ['foldx_Repair', 'foldx_Build', 'combine_results']
    # process_steps = ['foldx_Build', 'combine_results']
    # process_steps = []

    # foldx path
    if opsys=='Darwin':
        foldx_abspath = os.path.abspath('../../../yasara/foldx_2025/foldx_20251231_mac')
    elif opsys=='Linux':
        foldx_abspath = os.path.abspath('../../../yasara/foldx_2025/foldx_20251231_linux')
    elif opsys=='Windows':
        foldx_abspath = os.path.abspath('../../../yasara/foldx_2025/foldx_20251231_windows.exe')

    # get mutations
    # input is a list of positions to mutate
    res_mut_dict = {}
    if input_fname.find('.txt')>-1:
        with open(input_dir + input_fname) as f:
            mutations = [mut.replace('\n','') for mut in f.readlines()]
            # only WT positions specified, not mutations
            if mutations[0][-1].isdigit():
                wildtype_list = mutations.copy()
                # mutate to all possible residues, if not specified
                for wt in wildtype_list:
                    res_mut_dict[wt] = [wt + aa for aa in aaList if wt[0] != aa]
                mutations = [item for sublist in list(res_mut_dict.values()) for item in sublist]
            # both WT and MT specified
            elif mutations[0][-1].isalpha():
                wildtype_list = []
                for mut in mutations:
                    wt = mut[:-1]
                    if wt not in wildtype_list:
                        wildtype_list.append(wt)
                        res_mut_dict[wt] = []
                    res_mut_dict[wt].append(mut)
    # input is the protein sequence --> mutate every position
    elif input_fname.find('.fasta')>-1:
        seq_dict = read_fasta(input_dir + input_fname)
        seq = list(seq_dict.values())[0]
        wildtype_list = []
        for i,wt_aa in enumerate(seq):
            wt = wt_aa+str(i+1)
            res_mut_dict[wt] = [wt+aa for aa in aaList if wt!=aa]
            wildtype_list.append(wt)
        mutations = [item for sublist in list(res_mut_dict.values()) for item in sublist]
    print('# of mutations:', len(mutations))
    print('# of positions to mutate:', len(wildtype_list))

    # get mutation groups
    num_groups = int(np.ceil(len(wildtype_list) / num_res_per_group))
    print('# of mutation groups to process:', num_groups)
    mutations_group_list = []
    for grp_num in range(num_groups):
        res_grp = wildtype_list[
                  grp_num * num_res_per_group:min((grp_num + 1) * num_res_per_group, len(wildtype_list))]
        mutations_group_sublist = []
        for res in res_grp:
            mutations_group_sublist += res_mut_dict[res]
        mutations_group_list.append(mutations_group_sublist)
        print(mutations_group_sublist)

    ##############################
    ## RUN STABILITY PROCESSING ##
    ##############################

    # WITHOUT MULTIPROCESSING
    if run_multiprocessing is None:
        for i, mutation in enumerate(mutations):
            yasara_pid_i = run_stability_pipeline(mutation, input_pdb_fname, process_steps,
                                   pdb_dir, swapRes_dir, repair_dir, build_dir, log_fname, receptor_molname,
                                   remove_existing_dir=remove_existing_dir, foldx_abspath=foldx_abspath, multiprocessing_proc_num=None)
            if i==0:
                yasara_pid = yasara_pid_i
        # exit yasara after processing all mutations
        if yasara_pid is not None: exit_program(yasara_pid)

    else:
        num_cpu = cpu_count()
        print("Number of cpu : ", num_cpu)
        pool = Pool(processes=min(int(num_cpu/2),run_multiprocessing))
        for proc_num, mutations_group_sublist in enumerate(mutations_group_list):
            args = (mutations_group_sublist, input_pdb_fname, process_steps, pdb_dir, swapRes_dir, repair_dir, build_dir, log_fname, receptor_molname,
                    remove_existing_dir, foldx_abspath, proc_num)
            pool.apply_async(stability_pipeline_wrapper, args=args)
            print('Spawned process # ' + str(proc_num) + '. Mutations:' + ', '.join(mutations_group_sublist))
            time.sleep(5)
        pool.close()
        # for proc_num, mutation in enumerate(mutations):
        #     args = (mutation, input_pdb_fname, process_steps, pdb_dir, swapRes_dir, repair_dir, build_dir, log_fname,
        #             remove_existing_dir, foldx_abspath, proc_num)
        #     pool.apply_async(run_stability_pipeline, args=args)
        #     print('Spawned process # ' + str(proc_num) + '. Mutation:' + str(mutation))
        #     time.sleep(3)
        # pool.close()
        try:
            pool.join()
            print('Joined multiprocessing pool.')
        except:
            print('Unable to join pool. Moving on...')
        print('Completed all ' + str(proc_num+1) + ' processes.')

    # combine files spawned
    if 'combine_results' in process_steps:
        print('Combining stability calculation results for all mutants...')
        stability_cols = ['total energy', 'Backbone Hbond', 'Sidechain Hbond', 'Van der Waals', 'Electrostatics',
                          'Solvation Polar', 'Solvation Hydrophobic', 'Van der Waals clashes', 'entropy sidechain',
                          'entropy mainchain', 'sloop_entropy', 'mloop_entropy', 'cis_bond', 'torsional clash',
                          'backbone clash', 'helix dipole', 'water bridge', 'disulfide', 'electrostatic kon',
                          'partial covalent bonds', 'energy Ionisation', 'Entropy Complex']
        cols = ['Mutation', 'DDG'] + [col+'_MT' for col in stability_cols] + [col+'_WT' for col in stability_cols]

        mutations_list = []
        DDG_list = []
        wt_features_list = []
        mt_features_list = []

        for root, dirs, files in os.walk(build_dir):
            for name in files:
                if name.startswith("Raw"):
                    mut = name.split("_")[1]
                    mutations_list.append(mut)
                    with open(os.path.join(root, name)) as f:
                        # parse raw data
                        lines = f.read().splitlines()
                        header = lines[8].split('\t')
                        data = [l.split('\t') for l in lines[9:]]
                        # get mutant features
                        data_mt = np.array([row[1:] for row in data[:5]]).astype(float)
                        data_mt_avg = np.mean(data_mt,axis=0)
                        mt_features_list.append(list(data_mt_avg))
                        # get wildtype features
                        data_wt = np.array([row[1:] for row in data[5:]]).astype(float)
                        data_wt_avg = np.mean(data_wt, axis=0)
                        wt_features_list.append(list(data_wt_avg))
                        # get DDG
                        DDG = data_mt_avg[0] - data_wt_avg[0]
                        DDG_list.append(DDG)


        txt_all_list = [','.join(cols)]
        for mut, DDG, mt_features, wt_features in zip(mutations_list, DDG_list, mt_features_list, wt_features_list):
            txt_all_list.append(','.join(list(map(str, [mut, DDG]+mt_features+wt_features))))

        if os.path.exists(output_dir + output_fname + '.csv'):
            write_mode = 'a'
            txt_all_list = txt_all_list[1:]
        else:
            write_mode = 'w'

        # get text string to write
        txt_all = '\n'.join(txt_all_list)
        txt_all = txt_all.replace('\n\n', '\n').replace(',\n', '\n')

        # update or save file
        with open(output_dir + output_fname, write_mode) as f:
            f.write(txt_all)

if __name__ == "__main__":  # confirms that the code is under main function
    main()
