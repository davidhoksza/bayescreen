#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Compute RDKit descriptors for molecules/fragments.

Usage:
    python rdkit_descriptors.py
        -i {path to JSON with molecules, output of extract_fragments}
        -o {path to output CSV file}
        --fragments Use fragments else use molecules.
                    Default is to use molecules.

This file can be also imported as a python script. In such case please
use the extract_fragments method.
"""

import os
import argparse
import logging
import json
import rdkit
import rdkit.Chem
from rdkit.Chem import Descriptors

__author__ = 'Petr Škoda'
__license__ = 'X11'
__email__ = 'skoda@ksi.mff.cuni.cz'

# region Descriptors definition

_names = [
    'MolWt',
    'HeavyAtomMolWt',
    'ExactMolWt',
    'NumValenceElectrons',
    'NumRadicalElectrons',
    'MaxEStateIndex',
    'MinEStateIndex',
    'MaxAbsEStateIndex',
    'MinAbsEStateIndex',
    'BalabanJ',
    'BertzCT',
    'Chi0',
    'Chi0n',
    'Chi0v',
    'Chi1',
    'Chi1n',
    'Chi1v',
    'Chi2n',
    'Chi2v',
    'Chi3n',
    'Chi3v',
    'Chi4n',
    'Chi4v',
    'EState_VSA1',
    'EState_VSA10',
    'EState_VSA11',
    'EState_VSA2',
    'EState_VSA3',
    'EState_VSA4',
    'EState_VSA5',
    'EState_VSA6',
    'EState_VSA7',
    'EState_VSA8',
    'EState_VSA9',
    'FractionCSP3',
    'HallKierAlpha',
    'HeavyAtomCount',
    'Ipc',
    'Kappa1',
    'Kappa2',
    'Kappa3',
    'LabuteASA',
    'MolLogP',
    'MolMR',
    'NHOHCount',
    'NOCount',
    'NumAliphaticCarbocycles',
    'NumAliphaticHeterocycles',
    'NumAliphaticRings',
    'NumAromaticCarbocycles',
    'NumAromaticHeterocycles',
    'NumAromaticRings',
    'NumHAcceptors',
    'NumHDonors',
    'NumHeteroatoms',
    'NumRotatableBonds',
    'NumSaturatedCarbocycles',
    'NumSaturatedHeterocycles',
    'NumSaturatedRings',
    'PEOE_VSA1',
    'PEOE_VSA10',
    'PEOE_VSA11',
    'PEOE_VSA12',
    'PEOE_VSA13',
    'PEOE_VSA14',
    'PEOE_VSA2',
    'PEOE_VSA3',
    'PEOE_VSA4',
    'PEOE_VSA5',
    'PEOE_VSA6',
    'PEOE_VSA7',
    'PEOE_VSA8',
    'PEOE_VSA9',
    'RingCount',
    'SMR_VSA1',
    'SMR_VSA10',
    'SMR_VSA2',
    'SMR_VSA3',
    'SMR_VSA4',
    'SMR_VSA5',
    'SMR_VSA6',
    'SMR_VSA7',
    'SMR_VSA8',
    'SMR_VSA9',
    'SlogP_VSA1',
    'SlogP_VSA10',
    'SlogP_VSA11',
    'SlogP_VSA12',
    'SlogP_VSA2',
    'SlogP_VSA3',
    'SlogP_VSA4',
    'SlogP_VSA5',
    'SlogP_VSA6',
    'SlogP_VSA7',
    'SlogP_VSA8',
    'SlogP_VSA9',
    'TPSA',
    'VSA_EState1',
    'VSA_EState10',
    'VSA_EState2',
    'VSA_EState3',
    'VSA_EState4',
    'VSA_EState5',
    'VSA_EState6',
    'VSA_EState7',
    'VSA_EState8',
    'VSA_EState9',
    'fr_Al_COO',
    'fr_Al_OH',
    'fr_Al_OH_noTert',
    'fr_ArN',
    'fr_Ar_COO',
    'fr_Ar_N',
    'fr_Ar_NH',
    'fr_Ar_OH',
    'fr_COO',
    'fr_COO2',
    'fr_C_O',
    'fr_C_O_noCOO',
    'fr_C_S',
    'fr_HOCCN',
    'fr_Imine',
    'fr_NH0',
    'fr_NH1',
    'fr_NH2',
    'fr_N_O',
    'fr_Ndealkylation1',
    'fr_Ndealkylation2',
    'fr_Nhpyrrole',
    'fr_SH',
    'fr_aldehyde',
    'fr_alkyl_carbamate',
    'fr_alkyl_halide',
    'fr_allylic_oxid',
    'fr_amide',
    'fr_amidine',
    'fr_aniline',
    'fr_aryl_methyl',
    'fr_azide',
    'fr_azo',
    'fr_barbitur',
    'fr_benzene',
    'fr_benzodiazepine',
    'fr_bicyclic',
    'fr_diazo',
    'fr_dihydropyridine',
    'fr_epoxide',
    'fr_ester',
    'fr_ether',
    'fr_furan',
    'fr_guanido',
    'fr_halogen',
    'fr_hdrzine',
    'fr_hdrzone',
    'fr_imidazole',
    'fr_imide',
    'fr_isocyan',
    'fr_isothiocyan',
    'fr_ketone',
    'fr_ketone_Topliss',
    'fr_lactam',
    'fr_lactone',
    'fr_methoxy',
    'fr_morpholine',
    'fr_nitrile',
    'fr_nitro',
    'fr_nitro_arom',
    'fr_nitro_arom_nonortho',
    'fr_nitroso',
    'fr_oxazole',
    'fr_oxime',
    'fr_para_hydroxylation',
    'fr_phenol',
    'fr_phenol_noOrthoHbond',
    'fr_phos_acid',
    'fr_phos_ester',
    'fr_piperdine',
    'fr_piperzine',
    'fr_priamide',
    'fr_prisulfonamd',
    'fr_pyridine',
    'fr_quatN',
    'fr_sulfide',
    'fr_sulfonamd',
    'fr_sulfone',
    'fr_term_acetylene',
    'fr_tetrazole',
    'fr_thiazole',
    'fr_thiocyan',
    'fr_thiophene',
    'fr_unbrch_alkane',
    'fr_urea'
]

_functions = [
    Descriptors.MolWt,
    Descriptors.HeavyAtomMolWt,
    Descriptors.ExactMolWt,
    Descriptors.NumValenceElectrons,
    Descriptors.NumRadicalElectrons,
    # http://www.rdkit.org/docs/api/rdkit.Chem.Descriptors-module.html
    Descriptors.MaxEStateIndex,
    Descriptors.MinEStateIndex,
    Descriptors.MaxAbsEStateIndex,
    Descriptors.MinAbsEStateIndex,
    Descriptors.BalabanJ,
    Descriptors.BertzCT,
    Descriptors.Chi0,
    Descriptors.Chi0n,
    Descriptors.Chi0v,
    Descriptors.Chi1,
    Descriptors.Chi1n,
    Descriptors.Chi1v,
    Descriptors.Chi2n,
    Descriptors.Chi2v,
    Descriptors.Chi3n,
    Descriptors.Chi3v,
    Descriptors.Chi4n,
    Descriptors.Chi4v,
    Descriptors.EState_VSA1,
    Descriptors.EState_VSA10,
    Descriptors.EState_VSA11,
    Descriptors.EState_VSA2,
    Descriptors.EState_VSA3,
    Descriptors.EState_VSA4,
    Descriptors.EState_VSA5,
    Descriptors.EState_VSA6,
    Descriptors.EState_VSA7,
    Descriptors.EState_VSA8,
    Descriptors.EState_VSA9,
    Descriptors.FractionCSP3,
    Descriptors.HallKierAlpha,
    Descriptors.HeavyAtomCount,
    Descriptors.Ipc,
    Descriptors.Kappa1,
    Descriptors.Kappa2,
    Descriptors.Kappa3,
    Descriptors.LabuteASA,
    Descriptors.MolLogP,
    Descriptors.MolMR,
    Descriptors.NHOHCount,
    Descriptors.NOCount,
    Descriptors.NumAliphaticCarbocycles,
    Descriptors.NumAliphaticHeterocycles,
    Descriptors.NumAliphaticRings,
    Descriptors.NumAromaticCarbocycles,
    Descriptors.NumAromaticHeterocycles,
    Descriptors.NumAromaticRings,
    Descriptors.NumHAcceptors,
    Descriptors.NumHDonors,
    Descriptors.NumHeteroatoms,
    Descriptors.NumRotatableBonds,
    Descriptors.NumSaturatedCarbocycles,
    Descriptors.NumSaturatedHeterocycles,
    Descriptors.NumSaturatedRings,
    Descriptors.PEOE_VSA1,
    Descriptors.PEOE_VSA10,
    Descriptors.PEOE_VSA11,
    Descriptors.PEOE_VSA12,
    Descriptors.PEOE_VSA13,
    Descriptors.PEOE_VSA14,
    Descriptors.PEOE_VSA2,
    Descriptors.PEOE_VSA3,
    Descriptors.PEOE_VSA4,
    Descriptors.PEOE_VSA5,
    Descriptors.PEOE_VSA6,
    Descriptors.PEOE_VSA7,
    Descriptors.PEOE_VSA8,
    Descriptors.PEOE_VSA9,
    Descriptors.RingCount,
    Descriptors.SMR_VSA1,
    Descriptors.SMR_VSA10,
    Descriptors.SMR_VSA2,
    Descriptors.SMR_VSA3,
    Descriptors.SMR_VSA4,
    Descriptors.SMR_VSA5,
    Descriptors.SMR_VSA6,
    Descriptors.SMR_VSA7,
    Descriptors.SMR_VSA8,
    Descriptors.SMR_VSA9,
    Descriptors.SlogP_VSA1,
    Descriptors.SlogP_VSA10,
    Descriptors.SlogP_VSA11,
    Descriptors.SlogP_VSA12,
    Descriptors.SlogP_VSA2,
    Descriptors.SlogP_VSA3,
    Descriptors.SlogP_VSA4,
    Descriptors.SlogP_VSA5,
    Descriptors.SlogP_VSA6,
    Descriptors.SlogP_VSA7,
    Descriptors.SlogP_VSA8,
    Descriptors.SlogP_VSA9,
    Descriptors.TPSA,
    Descriptors.VSA_EState1,
    Descriptors.VSA_EState10,
    Descriptors.VSA_EState2,
    Descriptors.VSA_EState3,
    Descriptors.VSA_EState4,
    Descriptors.VSA_EState5,
    Descriptors.VSA_EState6,
    Descriptors.VSA_EState7,
    Descriptors.VSA_EState8,
    Descriptors.VSA_EState9,
    Descriptors.fr_Al_COO,
    Descriptors.fr_Al_OH,
    Descriptors.fr_Al_OH_noTert,
    Descriptors.fr_ArN,
    Descriptors.fr_Ar_COO,
    Descriptors.fr_Ar_N,
    Descriptors.fr_Ar_NH,
    Descriptors.fr_Ar_OH,
    Descriptors.fr_COO,
    Descriptors.fr_COO2,
    Descriptors.fr_C_O,
    Descriptors.fr_C_O_noCOO,
    Descriptors.fr_C_S,
    Descriptors.fr_HOCCN,
    Descriptors.fr_Imine,
    Descriptors.fr_NH0,
    Descriptors.fr_NH1,
    Descriptors.fr_NH2,
    Descriptors.fr_N_O,
    Descriptors.fr_Ndealkylation1,
    Descriptors.fr_Ndealkylation2,
    Descriptors.fr_Nhpyrrole,
    Descriptors.fr_SH,
    Descriptors.fr_aldehyde,
    Descriptors.fr_alkyl_carbamate,
    Descriptors.fr_alkyl_halide,
    Descriptors.fr_allylic_oxid,
    Descriptors.fr_amide,
    Descriptors.fr_amidine,
    Descriptors.fr_aniline,
    Descriptors.fr_aryl_methyl,
    Descriptors.fr_azide,
    Descriptors.fr_azo,
    Descriptors.fr_barbitur,
    Descriptors.fr_benzene,
    Descriptors.fr_benzodiazepine,
    Descriptors.fr_bicyclic,
    Descriptors.fr_diazo,
    Descriptors.fr_dihydropyridine,
    Descriptors.fr_epoxide,
    Descriptors.fr_ester,
    Descriptors.fr_ether,
    Descriptors.fr_furan,
    Descriptors.fr_guanido,
    Descriptors.fr_halogen,
    Descriptors.fr_hdrzine,
    Descriptors.fr_hdrzone,
    Descriptors.fr_imidazole,
    Descriptors.fr_imide,
    Descriptors.fr_isocyan,
    Descriptors.fr_isothiocyan,
    Descriptors.fr_ketone,
    Descriptors.fr_ketone_Topliss,
    Descriptors.fr_lactam,
    Descriptors.fr_lactone,
    Descriptors.fr_methoxy,
    Descriptors.fr_morpholine,
    Descriptors.fr_nitrile,
    Descriptors.fr_nitro,
    Descriptors.fr_nitro_arom,
    Descriptors.fr_nitro_arom_nonortho,
    Descriptors.fr_nitroso,
    Descriptors.fr_oxazole,
    Descriptors.fr_oxime,
    Descriptors.fr_para_hydroxylation,
    Descriptors.fr_phenol,
    Descriptors.fr_phenol_noOrthoHbond,
    Descriptors.fr_phos_acid,
    Descriptors.fr_phos_ester,
    Descriptors.fr_piperdine,
    Descriptors.fr_piperzine,
    Descriptors.fr_priamide,
    Descriptors.fr_prisulfonamd,
    Descriptors.fr_pyridine,
    Descriptors.fr_quatN,
    Descriptors.fr_sulfide,
    Descriptors.fr_sulfonamd,
    Descriptors.fr_sulfone,
    Descriptors.fr_term_acetylene,
    Descriptors.fr_tetrazole,
    Descriptors.fr_thiazole,
    Descriptors.fr_thiocyan,
    Descriptors.fr_thiophene,
    Descriptors.fr_unbrch_alkane,
    Descriptors.fr_urea
]


# endregion Descriptors definition

def create_parent_directory(path):
    """Create directory if it does not exists.

    :param path:
    :return:
    """
    dir_name = os.path.dirname(path)
    if not os.path.exists(dir_name) and not dir_name == "":
        os.makedirs(dir_name)


def _read_configuration():
    """Get and return application settings.

    :return:
    """
    parser = argparse.ArgumentParser(
        description='Compute RDKit descriptors for given'
                    'molecules/fragments.')
    parser.add_argument('-i', type=str, dest='input',
                        help='input JSON file',
                        required=True)
    parser.add_argument('-o', type=str, dest='output',
                        help='output CSV file', required=True)
    parser.add_argument('--fragments', dest='fragments',
                        help='use fragments instead of molecules',
                        action='store_true', required=False)

    return vars(parser.parse_args())


def compute_descriptors(input_file, output_file, use_fragments):
    """Compute descriptors for molecules/fragments in given input file.

    :param input_file:
    :param output_file:
    :param use_fragments: If true use fragments instead of molecules.
    :return: Summary object.
    """
    with open(input_file, 'r') as stream:
        data = json.load(stream)
    create_parent_directory(output_file)
    # Gather data.
    smiles_set = set()
    if use_fragments:
        for molecule in data:
            for fragment in molecule['fragments']:
                if not fragment['smiles'] in smiles_set:
                    smiles_set.add(fragment['smiles'])
    else:
        for molecule in data:
            if not molecule['smiles'] in smiles_set:
                smiles_set.add(molecule['smiles'])
    # Compute and write descriptors.
    sanitize_operation = rdkit.Chem.SanitizeFlags.SANITIZE_ALL ^ \
                         rdkit.Chem.SanitizeFlags.SANITIZE_KEKULIZE
    number_of_invalid = 0
    with open(output_file, 'w') as stream:
        stream.write('smiles,')
        stream.write(','.join(_names))
        stream.write('\n')
        counter = 0
        counter_step = int(len(smiles_set) / 10)
        for smiles in smiles_set:
            if counter % counter_step == 0:
                logging.info('%d/%d', counter, len(smiles_set))
            counter += 1
            # SMILES.
            stream.write('"')
            stream.write(smiles)
            stream.write('",')
            # Construct molecule, compute and write properties.
            molecule = rdkit.Chem.MolFromSmiles(str(smiles), sanitize=False)
            # Do not kekulize molecule.
            rdkit.Chem.SanitizeMol(molecule, sanitizeOps=sanitize_operation)
            #
            if molecule is None:
                logging.error('Invalid molecule detected: %s', smiles)
                number_of_invalid += 1
                continue
            stream.write(','.join([str(fnc(molecule)) for fnc in _functions]))
            stream.write('\n')
    # Log nad return summary.
    logging.info('Invalid molecules: %d/%d', number_of_invalid, len(smiles_set))
    return {
        'number_of_invalid': number_of_invalid,
        'total': len(smiles_set)
    }


def _main():
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s [%(levelname)s] %(module)s - %(message)s',
        datefmt='%H:%M:%S')
    configuration = _read_configuration()
    #
    use_fragments = 'fragments' in configuration and configuration['fragments']
    compute_descriptors(configuration['input'], configuration['output'],
                        use_fragments)


if __name__ == '__main__':
    _main()
