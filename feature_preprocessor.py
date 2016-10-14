#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Processing fragments' features CSV files. The CSV files are in the format
 returned by biochem-tools developed by Petr Skoda. Each line consists of a fragment SMILES
and a list of features. The header lists features names.
"""

import csv
import math
import numpy as np
import common
import logging

__author__ = "David Hoksza"
__email__ = "david.hoksza@mff.cuni.cz"
__license__ = 'X11'

# import matplotlib.pyplot as plt

def clusters_to_join(clusters, corr_matrix, corr_threshold):
    """
    Finds a pairs of clusters where every pair of compounds is correlated.

    :param clusters: List of existing clusters.
    :param corr_matrix: Correlation matrix storing correlation information for every pair of molecules.
    :param corr_threshold: The minimum correlation value for a pair of molecules to be considered correlated.
    :return: Pair (list) of correlated clusters or [-1,-1]
    """
    for i in range(len(clusters)):
        for j in range(i+1, len(clusters)):
            correlated = True
            for corr_i in clusters[i]:
                if not correlated: break
                for corr_j in clusters[j]:
                    if abs(corr_matrix[corr_i][corr_j]) < corr_threshold:
                        correlated = False
                        break
            if correlated: return [i, j]

    return [-1, -1]


def process(fn_csv_actives, fn_csv_inactives, log_file_name, corr_threshold):
    """
    The function takes input CSV files with fragments and their corresponding feature values, joins them
     into single list and:

    1. Removes features which are constant
    2. Imputes missing values by medians of the values for each feature.
    3. Computes correlation matrix for every pair of features
    4. Identifies clusters of correlated features.
    5. For each cluster, randomlly picks one feature to keep.

    :param fn_csv_actives: CSV file with fragments from active molecules and corresponding features.
    :param fn_csv_inactives: CSV file with fragments from inactive molecules and corresponding features.
    :param log_file_name: Log file name where to store information about correlation clusters.
    :param corr_threshold: The minimum correlation value for a pair of molecules to be considered correlated.
    :return: Dictionary with the remaining, preprocessed, features.
    """
    logging.info("Preprocessing descriptors (imputation, correlated descriptors removal)...")
    features = []
    feature_names = []
    fragment_names = []
    fragment_affiliations = {}
    #Reading features from CSVs
    for ix in range(0,2):
        flag = ""
        if ix == 0:
            fn = fn_csv_actives
            flag = "a"
        else:
            fn = fn_csv_inactives
            flag = "i"
        for row in csv.reader(open(fn, "r")):
            if not row: continue
            if len(features) == 0:
                for feature_name in row[1:]:
                    feature_names.append(feature_name)
                    features.append([])
            else:
                if row[0] not in fragment_names:
                    fragment_names.append(row[0])
                    fragment_affiliations[row[0]] = [flag]
                    for i in range(1, len(row)):
                        features[i-1].append(row[i])
                else:
                    fragment_affiliations[row[0]].append(flag)

    #Convert to numbers
    features = [[common.to_float(y) for y in x] for x in features]

    #Removal of constant features
    for i in range(len(features)-1, -1, -1):
        if np.isnan(features[i]).all() or np.nanmax(features[i]) - np.nanmin(features[i]) == 0:
            del features[i]
            del feature_names[i]

    #Get medians and impute them
    imputation_values = []
    for i in range(len(features)):
        median = np.median(np.asarray([x for x in features[i] if not math.isnan(x)]))
        imputation_values.append(median)
        features[i] = [median if math.isnan(x) else x for x in features[i]]

    corr_matrix = [[0 for x in range(len(features))] for x in range(len(features))]
    for i in range(len(features)):
        a = np.asarray(features[i])
        for j in range(i, len(features)):
            c = np.corrcoef(a, np.asarray(features[j]))[0, 1]
            corr_matrix[i][j] = c
            corr_matrix[j][i] = c

    #Add each feature into its own cluster
    clusters = []
    for i in range(len(features)):
        clusters.append([i])

    while True:
        [i, j] = clusters_to_join(clusters, corr_matrix, corr_threshold)
        if i == -1: break
        clusters[i] += clusters[j]
        clusters.pop(j)

    #Print clusters
    if log_file_name is not None and log_file_name != "":
        with open(log_file_name, "w") as fl:
            for cluster in clusters:
                str_cluster = feature_names[cluster[0]] + ": "
                for ix_feature in cluster:
                    str_cluster += " " + feature_names[ix_feature]
                fl.write(str_cluster + "\n")

    #Leave out all but first features of each cluster from the feature matrix
    ixs_features_to_remove = sorted([j for i in clusters for j in i[1:]])
    for ix in reversed(ixs_features_to_remove):
        feature_names.pop(ix)
        features.pop(ix)

    #Print out the resulting matrix
    # with open(out_file_name, "w") as fo:
    #     line = "Name"
    #     for name in feature_names: line += ",{}".format(name)
    #     fo.write(line + "\n")
    #     for i in range(len(features[0])):
    #         line = fragment_names[i]
    #         for feature in features: line += ",{}".format(feature[i])
    #         fo.write(line + "\n")

    return {"features_names": feature_names,
            "fragments_names": fragment_names,
            "features_vectors": features,
            "affiliations": fragment_affiliations,
            "imputation_values": imputation_values}

    #fig, ax = plt.subplots()
    #heatmap = ax.pcolor(corrMatrix, cmap=plt.cm.Blues, vmin=0, vmax=2)

    #ax.set_xticklabels(feature_names, minor=False)
    #ax.set_yticklabels(feature_names, minor=False)
    #plt.show()
