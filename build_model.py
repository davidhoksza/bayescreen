#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

The purpose of this script is to compute the probabilities comprising the Bayes model of activity.
The steps carried out to build the model from known active and inactive molecules, passed on input, include:

* Extract fragmetns from the input molecules (main function)
* Generate descriptors for all the fragments (main function)
* Preprocess the features (removal of correlated features, normalization, binning) (main adn build_model function)
* Computation of the probabilities which comprise the model (build_model function)
* Writing out the model in JSON format (main function)

"""

import argparse
import numpy as np
import csv
import json
import sys
import math
from rdkit import Chem
import common
import logging

import biochem_tools
from feature_preprocessor import process as preprocess_features

__author__ = "David Hoksza"
__email__ = "david.hoksza@mff.cuni.cz"
__license__ = 'X11'


def build_model(fn_actives_json, fn_inactives_json, features, cnt_bins=100):
    """
    Builds an activity model from active and inactive feature vectors. These feature vectors are expected
    to be non-correlated.

    :param fn_actives_json: JSON with the active fragments; for each molecules there is a list of its fragments.
    :param fn_inactives_json: JSON with the inactive fragments; for each molecules there is a list of its fragments.
    :param features: Dictionary with (fragment; feature vector pairs)
    :param cnt_bins: Number of bins to be used for binning
    :return: JSON file with the model
    """
    logging.info("Building model...")
    feature_matrix = np.array(features["features_vectors"]).transpose()

    #Normalize
    mins = np.nanmin(feature_matrix, axis=0)
    maxs = np.nanmax(feature_matrix, axis=0)

    feature_matrix = (feature_matrix - mins) / (maxs-mins)

    #Serialize the normalized features vectors matrix
    # with open("normalized.features.csv", "w") as f:
    #     line = "Name"
    #     for name in features["features_names"]: line += ",{}".format(name)
    #     f.write(line + "\n")
    #
    #     line = "mins"
    #     for m in mins: line += ",{}".format(m)
    #     f.write(line + "\n")
    #
    #     line = "maxs"
    #     for m in maxs: line += ",{}".format(m)
    #     f.write(line + "\n")
    #
    #     for ix in range(len(features["fragments_names"])):
    #         line = features["fragments_names"][ix]
    #         for feature in feature_matrix[ix]: line += ",{}".format(feature)
    #         f.write(line + "\n")


    train_features = {"actives": [],
                      "inactives": []}  # stores feature vectors irrespective of which molecule they come from

    fns = [fn_actives_json, fn_inactives_json]
    activity = ["actives", "inactives"]
    for ix in range(0, 2):
        with open(fns[ix], "r") as f:
            molecules = json.load(f)
            for mol in molecules:
                processed_fragments = []
                for fragment in [frag["smiles"] for frag in mol["fragments"]]:
                    if fragment not in processed_fragments:
                        processed_fragments.append(fragment)
                        features_vector = feature_matrix[features["fragments_names"].index(fragment)].tolist()
                        train_features[activity[ix]].append(features_vector)

    #Obtain histograms (cnt_bins) for every feature in actives and inactives
    histograms = {"actives": [], "inactives": []}
    aux_features = np.array(train_features["actives"])
    for ix in range(0, len(features["features_names"])):
        histograms["actives"].append(np.histogram(aux_features[:, ix], bins=cnt_bins, range=(0, 1))[0])
    aux_features = np.array(train_features["inactives"])
    for ix in range(0, len(features["features_names"])):
        histograms["inactives"].append(np.histogram(aux_features[:, ix], bins=cnt_bins, range=(0, 1))[0])

    #Compute the probabilities
    probs = {}
    probs["feature_value_in_actives"] = [] #P(f_i=X|A)
    for ix in range(0, len(features["features_names"])):
    # for ix in significant_features:
        probs["feature_value_in_actives"].append(((histograms["actives"][ix]+1)/float(len(train_features["actives"])+1)).tolist())
    probs["feature_value_in_inactives"] = [] #P(f_i=X|I)
    for ix in range(0, len(features["features_names"])):
    # for ix in significant_features:
        probs["feature_value_in_inactives"].append(((histograms["inactives"][ix]+1)/float(len(train_features["inactives"])+1)).tolist())
    probs["feature_value"] = [] #P(f_i=X)
    for ix in range(0, len(features["features_names"])):
    # for ix in significant_features:
        probs["feature_value"].append(((histograms["actives"][ix] + histograms["inactives"][ix] + 1) / float((len(train_features["actives"]) + len(train_features["inactives"]) + 1))).tolist())
    probs["active"] = len(train_features["actives"])/float(len(train_features["actives"]) + len(train_features["inactives"])) #P(A)
    probs["inactive"] = len(train_features["inactives"])/float(len(train_features["actives"]) + len(train_features["inactives"])) #P(I)

    model = {
        "probabilities": probs,
        "cnt_bins": cnt_bins,
        "features_names": features["features_names"],
        "fragment_types": args.fragment_types,
        "features_generator": args.descriptors_generator,
        "path_to_padel": args.padel_root,
        "normalization": {
            "mins": mins.tolist(),
            "maxs": maxs.tolist(),
            "imputation_values": features["imputation_values"]
        }
    }

    return json.dumps(model, indent=4)


def main():
    common.init_logging()
    [fn_actives_json, fn_inactives_json] = common.fragments_extraction([args.actives, args.inactives], args.fragment_types)
    [fn_actives_csv, fn_inactives_csv] = common.descriptors_extraction([fn_actives_json, fn_inactives_json],
                                                                       args.descriptors_generator,
                                                                       args.padel_root)
    features = preprocess_features(fn_actives_csv, fn_inactives_csv, args.log, args.corr_threshold)
    model = build_model(fn_actives_json, fn_inactives_json, features, args.cnt_bins)

    with common.open_file(args.model, "w") as f:
        f.write(model)

    if args.clean:
        common.delete_files([fn_actives_json, fn_inactives_json, fn_actives_csv, fn_inactives_csv])


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--actives",
                        required=True,
                        metavar='FILE',
                        type=str,
                        help="Known active molecules in SDF or SMILES format")
    parser.add_argument("-i", "--inactives",
                        required=True,
                        metavar = 'FILE',
                        type=str,
                        help="Known inactive molecules in SDF or SMILES format")
    parser.add_argument("-m", "--model",
                        required=True,
                        metavar='FILE',
                        type=str,
                        help="Output JSON file containing the resulting model.")
    parser.add_argument("-t", "--file-type",
                        help="Optional type of input files ('sdf', 'smi'). By default, the type is determined"
                             "from the filename extension, but can be overriden using this argument. If the"
                             "file has no extension, sdf will be used as default.",
                        default="sdf")
    parser.add_argument("-f", "--fragment-types",
                        help="Optional comma separated list of fragment types to extract (eg. \"tt.3,ecfp.2\")",
                        default="ecfp.2")
    parser.add_argument("-g", "--descriptors-generator",
                        help="Generator to be used to obtain fragments descriptors. Allowed values are [rdkit|padel]",
                        choices=['rdkit', 'padel'],
                        default="rdkit")
    parser.add_argument("-r", "--padel-root",
                        help="Path to PaDEL descriptors generator. Required when -g is set to 'padel'.")
    parser.add_argument("-l", "--log",
                        help="Log file name")
    parser.add_argument("-c", "--clean",
                        action='store_true',
                        help="Delete the fragment and features files.")
    parser.add_argument("--corr-threshold",
                        default=0.7,
                        type=float,
                        help="Features with correlation above this threshold are considered correlated.")
    parser.add_argument("--cnt-bins",
                        default="100",
                        type=int,
                        help="Number of bins when binning feature values")

    args = parser.parse_args()

    if args.descriptors_generator == "padel" and (args.padel_root is None ):
        raise argparse.ArgumentTypeError("--padel-root not set.")
    if args.descriptors_generator == "padel" and (args.padel_root is None ):
        raise argparse.ArgumentTypeError("--padel-root not set.")

    main()
