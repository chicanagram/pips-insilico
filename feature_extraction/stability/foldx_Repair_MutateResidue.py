# -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 15:11:37 2021

@author: Jhoann
"""

import subprocess
import os, sys

mapping = {
    'A': 'Ala',
    'H': 'His',
    'Y': 'Tyr',
    'R': 'Arg',
    'T': 'Thr',
    'K': 'Lys',
    'M': 'Met',
    'D': 'Asp',
    'N': 'Asn',
    'C': 'Cys',
    'Q': 'Gln',
    'E': 'Glu',
    'G': 'Gly',
    'I': 'Ile',
    'L': 'Leu',
    'F': 'Phe',
    'P': 'Pro',
    'S': 'Ser',
    'W': 'Trp',
    'V': 'Val'
    }
def get_mapping(mapping, res):
    # Map residue initial to code
    result = mapping[res]
    return result

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

# perform foldx RepairPDB
def foldxRepair(pdbs, output_dir, input_dir, foldx_path):
    print('Started Foldx Repair PDB')
    indv_list_fname = os.path.join(output_dir, 'pdbList.txt')
    with open(indv_list_fname, 'w') as output:
        if isinstance(pdbs, list):
            output.write('{0};\n'.format(';\n'.join(pdbs)))
        else:
            output.write('{0};\n'.format(pdbs))

    cmd = [
        foldx_path,
        '--command=RepairPDB',
        '--pdb-list={0}'.format(indv_list_fname),
        '--pdb-dir={0}'.format(input_dir),
        '--out-pdb=true',
        '--output-dir={0}'.format(output_dir)]

    foldx_args = " ".join(cmd)
    print(foldx_args)
    # subprocess.run(foldx_args, shell=True, encoding="utf8")
    subprocess.call(cmd)
    print('Ended Foldx Repair PDB')


def foldxBuild(mutation, pdb, output_dir, input_dir, foldx_path, receptor_molname='R'):
    reverseMut = []
    # convert single mutation to list of that mutation
    if not isinstance(mutation, list):
        mutation = [mutation]
    for mutgrp in mutation:
        WT, pos, MT = split_mutation(mutgrp, aa_letter_representation=True)
        revMut = MT + receptor_molname + str(pos) + WT
        reverseMut.append(revMut)
    print(reverseMut)

    print('Started Foldx BuildModel')
    indv_list_fname = os.path.join(output_dir, 'individual_list.txt')
    with open(indv_list_fname, 'w') as output:
        if isinstance(reverseMut, list):
            output.write('{0};\n'.format(','.join(reverseMut)))
        else:
            output.write('{0};\n'.format(reverseMut))

    cmd = [
            foldx_path,
           '--command=BuildModel',
           '--pdb={0}_Repair.pdb'.format(pdb),
           '--mutant-file={0}'.format(indv_list_fname),
           '--pdb-dir={0}'.format(input_dir),
           '--ionStrength=0.05',
           '--pH=7',
           '--vdwDesign=2', 
           '--out-pdb=true',
           '--numberOfRuns=5', 
           '--temperature=298',
           '--moveNeighbours=true', 
           '--output-dir={0}'.format(output_dir)]
    
    foldx_args = " ".join(cmd)
    # subprocess.run(foldx_args, shell=True, encoding="utf8")
    subprocess.call(cmd)
    print('Ended Foldx BuildModel')