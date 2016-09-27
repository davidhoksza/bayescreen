#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Extract fragments for molecules in given SDF files and save then into JSON.

Path and circular fragments use different indexing method, so they may collide.

Usage:
    python extract_fragments.py
        -i {input file or directory with input files}
        -o {path to output}
        -f {optional, comma separated list of fragment types to extract}

        -t {type of input files, 'sdf', 'smi'. Default is 'sdf'}
        --kekule {generated kekule form of SMILES for fragments}
        --isomeric {put stereochemistry information into fragments SMILES}

Fragments type:
    - tt.{SIZE}
    - ecfp.{SIZE}
where {SIZE} should be replaced by required fragment size. Usage example:
    tt.3,ecfp.2
default value:
    tt.3

Kekule smiles form has no aromatic bonds. Use of --kekule option thus may
reduce the number of generated unique fragments.

This file can be also imported as a python script. In such case please
use the extract_fragments method.
"""

import os
import argparse
import logging
import json
import rdkit
import rdkit.Chem
from rdkit.Chem import AllChem
import rdkit.Chem.AtomPairs.Utils

__author__ = 'Petr Škoda'
__license__ = 'X11'
__email__ = 'skoda@ksi.mff.cuni.cz'

# region Path fragments

atom_code = {
    'bits': {
        'type': 4,
        'pi': 2,
        'branch': 4,
        'total': 10
    }
}


def get_atom_code(atom, branch_subtract):
    # Constants;
    num_type_bits = atom_code['bits']['type']
    num_pi_bits = atom_code['bits']['pi']
    num_branch_bits = atom_code['bits']['branch']

    # code = typeIdx | numPiElectrons | numBranches

    max_num_branches = (1 << num_branch_bits) - 1
    max_num_pi = (1 << num_pi_bits) - 1
    # Original publication use :
    #                   [5, 6, 7, 8, 9, 14, 15, 16, 17, 33, 34, 35, 53]
    # RDKit use:
    # We must add trailing zero as we need 16 elements in the array
    # for atom_code.bits.type equal 4.
    atom_number_types = [5, 6, 7, 8, 9, 14, 15, 16, 17, 33, 34, 35, 51, 52, 43,
                         0]
    # Number of non-hydrogen? neighbor
    if atom.GetDegree() > branch_subtract:
        num_branches = atom.GetDegree() - branch_subtract
    else:
        num_branches = 0
    code = num_branches % max_num_branches

    # Number of bonding pi-electrons.
    n_pi = rdkit.Chem.AtomPairs.Utils.NumPiElectrons(atom) % max_num_pi
    code |= n_pi << num_branch_bits

    # If atom.getAtomicNum() is in atomNumberTypes then return
    # exact match. Otherwise return smallest bigger value.
    type_idx = 0
    n_types = 1 << num_type_bits;
    while type_idx < n_types:
        if atom_number_types[type_idx] == atom.GetAtomicNum():
            break
        elif atom_number_types[type_idx] > atom.GetAtomicNum():
            type_idx = n_types
            break
        else:
            type_idx += 1

    # Make sure we do not point outside the array.
    if type_idx == n_types:
        type_idx -= 1

    # Atom type.
    code |= type_idx << (num_branch_bits + num_pi_bits);

    return code


def score_path(molecule, path, size):
    codes = [None] * size
    for i in range(size):
        if i == 0 or i == (size - 1):
            sub = 1
        else:
            sub = 2
        # We use this branch airways as we do not use custom atomCodes.
        codes[i] = get_atom_code(molecule.GetAtomWithIdx(path[i]), sub)

    # We scan the vector for both sides, we want to make sure that
    # the begging is less or equal to the end.

    # "canonize" the code vector:
    beg = 0
    end = len(codes) - 1
    while beg < end:
        if codes[beg] > codes[end]:
            codes.reverse()
            break
        elif codes[beg] == codes[end]:
            beg += 1
            end -= 1
        else:
            break

    # Just add all together.
    accum = 0
    for i in range(size):
        accum |= (codes[i]) << (atom_code['bits']['total'] * i)
    return accum


def extract_path_fragments(molecule, size, options):
    output = []
    pattern = rdkit.Chem.MolFromSmarts('*' + ('~*' * (size - 1)))
    for atoms in molecule.GetSubstructMatches(pattern):
        smiles = rdkit.Chem.MolFragmentToSmiles(
            molecule, atomsToUse=list(atoms),
            kekuleSmiles=options['kekule'],
            isomericSmiles=options['isomeric'])
        output.append({
            'smiles': smiles,
            'index': score_path(molecule, atoms, size),
            'type': 'TT',
            'size': size
        })
    return output


# endregion

# region Circular fragments

def extract_neighbourhood_fragments(molecule, size, options):
    """Extract and return circular fragments.

    :param molecule:
    :param size:
    :param options:
    :return:
    """
    output = []
    info = {}
    AllChem.GetMorganFingerprint(molecule, radius=size, bitInfo=info)
    for element in info:
        for item in info[element]:
            # item = [rooted atom, radius]
            if item[1] < size:
                continue
            # assemble fragments into atom
            env = rdkit.Chem.FindAtomEnvironmentOfRadiusN(
                molecule, item[1], item[0])
            atoms = set()
            for bidx in env:
                atoms.add(molecule.GetBondWithIdx(bidx).GetBeginAtomIdx())
                atoms.add(molecule.GetBondWithIdx(bidx).GetEndAtomIdx())
            # check if we have some atoms
            if len(atoms) > 0:
                try:
                    # kekuleSmiles - we may lost some information
                    # about aromatic atoms, but if we do not kekulize
                    # we can get invalid smiles
                    smiles = rdkit.Chem.MolFragmentToSmiles(
                        molecule, atomsToUse=list(atoms), bondsToUse=env,
                        rootedAtAtom=item[0], kekuleSmiles=options['kekule'],
                        isomericSmiles=options['isomeric'])
                except Exception:
                    logging.exception('Invalid fragment detected.')
                    logging.info('Molecule: %s', molecule.GetProp('_Name'))
                    logging.info('Atoms: %s', ','.join([str(x) for x in atoms]))
                output.append({
                    'smiles': smiles,
                    'index': element,
                    'type': 'ECFP',
                    'size': size
                })
    return output


# endregion

def extract_fragments_from_molecule(molecule, types, options):
    """Return fragments for given molecule.

    :param molecule:
    :param types: Types of fragments to extract.
    :param options
    :return:
    """
    output = []
    for item in types:
        if item['name'] == 'tt':
            output.extend(extract_path_fragments(
                molecule, item['size'], options))
        elif item['name'] == 'ecfp':
            output.extend(extract_neighbourhood_fragments(
                molecule, item['size'], options))
    return output


def _read_configuration():
    """Get and return application settings.

    :return:
    """
    parser = argparse.ArgumentParser(
        description='Extract molecular fragments. '
                    'See file header for more details.')
    parser.add_argument('-i', type=str, dest='input', required=True)
    parser.add_argument('-o', type=str, dest='output', required=True)
    parser.add_argument('-f', type=str, dest='fragments', required=False)
    parser.add_argument('-t', type=str, dest='input_type', default='sdf')
    parser.add_argument('--recursive', dest='recursive', action='store_true',
                        required=False)
    parser.add_argument('--kekule', dest='kekule',
                        action='store_true', required=False)
    parser.add_argument('--isomeric', dest='isomeric',
                        action='store_true', required=False)

    configuration = vars(parser.parse_args());

    if 'fragments' not in configuration or configuration['fragments'] is None:
        configuration['fragments'] = 'tt.3'

    # Parse fragment types.
    parsed_types = []
    for item in configuration['fragments'].split(','):
        item_split = item.split('.')
        if not len(item_split) == 2:
            logging.error('Invalid fragment type: %s', item)
            logging.info('  Expected format {TYPE}.{SIZE}')
            exit(1)
        parsed_types.append({
            'name': item_split[0],
            'size': int(item_split[1])
        })
    configuration['fragments'] = parsed_types
    configuration['input_type'] = configuration['input_type'].lower()
    return configuration


def load_sdf(path):
    """Generate molecules from SDF file.

    :param path:
    :param types:
    """
    logging.info('Loading (SDF): %s' % path)
    for molecule in rdkit.Chem.SDMolSupplier(path):
        if molecule is None:
            logging.error('Invalid molecule detected.')
            continue
        yield molecule


def load_smi(path):
    """Generate molecules from SMI file.

    :param path:
    :return:
    """
    logging.info('Loading (SMI): %s' % path)
    with open(path, 'r') as stream:
        for line in stream:
            line = line.strip()
            molecule = rdkit.Chem.MolFromSmiles(line)
            if molecule is None:
                logging.error('Invalid molecule detected.')
                continue
            # Molecules created from SMILES does not have any name,
            # so we use the SMILES as a name.
            molecule.SetProp('_Name', line)
            yield molecule


def recursive_scan_for_input(path, recursive, extension):
    """Perform recursive scan for input files.

    :param path:
    :param recursive
    :param extension
    :return:
    """
    result = []
    for file_name in os.listdir(path):
        file_path = path + '/' + file_name
        if os.path.isdir(file_path):
            if recursive:
                result.extend(recursive_scan_for_input(
                    file_path, recursive, extension))
        elif os.path.isfile(file_path) \
                and file_name.lower().endswith(extension):
            result.append(file_path)
    return result


def append_object_to_json(output_stream, item, holder):
    """Write given molecule as a JSON into stream.

    Optionally put separator before the record based on 'holder'.
    :param output_stream:
    :param item: Item to append to JSON file.
    :param holder: Object shared by all calls of this method on the same stream.
    :return:
    """
    if holder['first']:
        holder['first'] = False
    else:
        output_stream.write(',')
    json.dump(item, output_stream)


def create_parent_directory(path):
    """Create directory if it does not exists.

    :param path:
    :return:
    """
    dir_name = os.path.dirname(path)
    if not os.path.exists(dir_name) and not dir_name == "":
        os.makedirs(dir_name)


_load_functions = {
    'sdf': load_sdf,
    'smi': load_smi
}


def extract_fragments(input_files, input_type, output_file, extraction_options):
    """Extract fragments from molecules and write them to output JSON file.

    The extraction_options['fragments'] must be a list with objects describing
    fragments to extract, see _read_configuration for more details.

    :param input_files: List of files with molecules.
    :param input_type: Type of input see _load_functions property.
    :param output_file: Path to output JSON file.
    :param extraction_options: See usage in _main for more information.
    :return: Object with summary about computation.
    """
    # The write_molecule_json need some static info.
    holder = {'first': True}
    # Count some statistics.
    total_fragments = 0
    #
    create_parent_directory(output_file)
    with open(output_file, 'w') as output_stream:
        output_stream.write('[')
        for path in input_files:
            for molecule in _load_functions[input_type](path):
                item = {
                    'name': molecule.GetProp('_Name'),
                    'smiles': rdkit.Chem.MolToSmiles(molecule),
                    'fragments': extract_fragments_from_molecule(
                        molecule, extraction_options['fragments'],
                        extraction_options)
                }
                total_fragments += len(item['fragments'])
                # Append to output.
                append_object_to_json(output_stream, item, holder)
        output_stream.write(']')
    # Log nad return summary.
    logging.info('Report')
    logging.info('\tfragments total: %d', total_fragments)
    return {
        'total_fragments': total_fragments
    }


def _main():
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s [%(levelname)s] %(module)s - %(message)s',
        datefmt='%H:%M:%S')
    configuration = _read_configuration()
    # Read files to load.
    if os.path.isdir(configuration['input']):
        input_files = recursive_scan_for_input(configuration['input'],
                                               configuration['recursive'],
                                               configuration['input_type'])
    else:
        input_files = [configuration['input']]
    # Prepare configuration for the extraction.
    extraction_options = {
        'kekule': configuration['kekule'],
        'isomeric': configuration['isomeric'],
        'fragments': configuration['fragments']
    }
    #
    extract_fragments(input_files, configuration['input_type'],
                      configuration['output'], extraction_options)


if __name__ == '__main__':
    _main()
