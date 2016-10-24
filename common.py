#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Analyzes Bayes activity model.

"""

import sys
import logging
import gzip
import numpy as np
import biochem_tools
import os
import fnmatch

__author__ = "David Hoksza"
__email__ = "david.hoksza@mff.cuni.cz"
__license__ = 'X11'


def find_files_recursively(directory, pattern):
    matches = []
    for root, dirnames, filenames in os.walk(directory):
        for filename in fnmatch.filter(filenames, pattern):
            matches.append(os.path.join(root, filename))
    return matches


def open_file(file_name, mode="r"):
    access_type = mode
    if sys.version_info >= (3,): access_type = mode + "t"
    if file_name.endswith("gz"):
        return gzip.open(file_name, access_type)
    else:
        return open(file_name, access_type)


def init_logging():
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s [%(levelname)s] %(module)s - %(message)s',
        datefmt='%H:%M:%S')


def to_float(x):
    try:
        a = float(x)
        if np.isinf(a): a = float('nan')
    except ValueError:
        return float('nan')
    else:
        return a


def fragments_extraction(ds_file_names, fragment_types):
    extraction_options = {
        'kekule': False,
        'isomeric': False,
        'fragments': fragment_types
    }

    parsed_types = []
    for item in fragment_types.split(','):
        item_split = item.split('.')
        if not len(item_split) == 2:
            logging.error('Invalid fragment type: %s', item)
            logging.info('  Expected format {TYPE}.{SIZE}')
            exit(1)
        parsed_types.append({
            'name': item_split[0],
            'size': int(item_split[1])
        })
    extraction_options['fragments'] = parsed_types

    fn_json = []
    for fn in ds_file_names:
        file_type = "sdf"
        if fn.endswith(".smi"):
            file_type = "smi"

        if fn.endswith(".smi") or fn.endswith(".sdf"):
            fn_json.append(fn[:-4] + ".frags.json")
        else:
            fn_json.append(fn + ".json")
        biochem_tools.extract_fragments([fn], file_type, fn_json[-1], extraction_options)

    return fn_json


def descriptors_extraction(json_file_names, descriptors_generator, padel_path):
    fn_csv = []
    for fn in json_file_names:
        fn_csv.append(fn[:-5] + ".csv")
        if descriptors_generator == "rdkit":
            biochem_tools.rdkit_compute_descriptors(fn, fn_csv[-1], True)
        elif descriptors_generator == "padel":
            biochem_tools.padel_compute_descriptors(fn, fn_csv[-1], True, padel_path)

    return fn_csv


def delete_files(file_names):
    for fn in file_names:
        os.remove(fn)
