#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Compute PaDEL descriptors for molecules/fragments.

Usage:
    python padel_descriptors.py
        -i {path to JSON with molecules, output of extract_fragments}
        -o {path to output csv file}
        -p {path to the PaDEL directory that contains PaDEL-Descriptor.jar}
        -f Compute for fragments else for molecules.


This file can also be used as a python script for import, in such case
please use the compute_descriptors method.
"""

import os
import argparse
import logging
import json
import subprocess

__author__ = 'Petr Škoda'
__license__ = 'X11'
__email__ = 'skoda@ksi.mff.cuni.cz'


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
        description='Compute PaDEL descriptors for given'
                    'molecules/fragments.')
    parser.add_argument('-i', type=str, dest='input',
                        help='input JSON file',
                        required=True)
    parser.add_argument('-o', type=str, dest='output',
                        help='output CSV file', required=True)
    parser.add_argument('-p', type=str, dest='padel',
                        help='PaDEL directory', required=True)
    parser.add_argument('-f', dest='fragments',
                        help='use fragments instead of molecules',
                        action='store_true', required=False)

    return vars(parser.parse_args())


def compute_descriptors(input_file, output_file, use_fragments, padel_path):
    """Compute descriptors for molecules/fragments in given input file.

    :param input_file:
    :param output_file:
    :param use_fragments: If true use fragments instead of molecules.
    :param padel_path: Path to PaDel.
    :return: Nothing.
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
    # Prepare data for PaDEL.
    padel_input = os.path.dirname(output_file) + '/PaDEL-temp.smi'
    with open(padel_input, 'w') as stream:
        for smiles in smiles_set:
            stream.write(smiles)
            stream.write('\t')
            stream.write(smiles)
            stream.write('\n')
    # Execute PaDEL.
    logging.info('Executing PaDEL ...')
    thread = subprocess.Popen(
        ['java', '-jar', '-Xmx1024m',
         padel_path + '/PaDEL-Descriptor.jar',
         # '-maxruntime', '5000',
         '-threads', '1',
         '-2d',
         '-dir', padel_input,
         '-file', output_file],
        shell=True)
    thread.wait()
    logging.info('Executing PaDEL ... done')
    os.remove(padel_input)
    # Return summary.
    return {
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
                        use_fragments, configuration['padel'])


if __name__ == '__main__':
    _main()
