#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Computes performance of a simulated virtual screening campaign where actives and inactives
are known in advance.

The performance characteristics are:
    AUC - Arean Under the ROC (Receiver Operating Characteristics) Curve
    EFn - Enrichment factor computed on top n% of the database
"""

import common
import argparse
import numpy as np
import os
from rdkit.ML.Scoring import Scoring


__author__ = "David Hoksza"
__email__ = "david.hoksza@mff.cuni.cz"
__license__ = 'X11'

EF_FRACTIONS = [0.005, 0.01, 0.02, 0.05]


def evaluate_pair(fn_actives, fn_inactives):
    """
    For a pair of files with actives and inactives scores computes
    AUC (area under the ROC curve) and EF (enrichment factor).
    :param fn_actives: Name of the file with the active molecules scores.
    :param fn_inactives: Name of the file with the inactive molecules scores.
    :return: A dictionary with "auc" and "ef" keys. "auc" contains a single value, "ef" contains
    a list of EFs corresponding to EF_FRACTIONS fractions.
    """
    ranking = []
    iteration = [{"fn": fn_actives, "activity": 1}, {"fn": fn_inactives, "activity": 0}]
    for it in iteration:
        with common.open_file(it["fn"]) as f:
            for line in f:
                s_line = line.split(":")
                ranking.append({"mol": s_line[0].strip(),
                                "score": common.to_float(s_line[1].strip()),
                                "activity": it["activity"]})
    ranking = sorted(ranking, key=lambda x: x["score"], reverse=True)
    return {
        "auc": Scoring.CalcAUC(ranking, 'activity'),
        "ef": Scoring.CalcEnrichment(ranking, 'activity', EF_FRACTIONS)
        }


def evaluate_directory(dir):
    """
    Evaluates multiple simulated screens. The file with actives score is expected to have "actives.out" suffix
     while the file with inactives scores is expected to have "inactives.out" suffix. The prefix needs to be
     identical for results from one campaign (e.g., CDK2_actives.out and CDK2_inactives.out)

    :param dir: Directory to be recursively searched for the actives and inactives scoring files.
    :return:
    """
    agg_results = {}
    for fn_actives in common.find_files_recursively(dir, "*actives.out"):
        fn_inactives = fn_actives.replace("actives", "inactives")
        performance = evaluate_pair(fn_actives, fn_inactives)
        print("AUC: {}\n{}".format(fn_actives.replace("actives", ""), performance["auc"]))
        for i in range(len(EF_FRACTIONS)):
            print("EF{}: {}".format(EF_FRACTIONS[i], performance["ef"][i]))

        s_fn = fn_inactives.replace("\\", "/").split("/")
        group = s_fn[args.group_index]
        target = s_fn[args.target_index]
        if group not in agg_results:
            agg_results[group] = {}
        if target not in agg_results[group]:
            agg_results[group][target] = []
        agg_results[group][target].append(performance)

    print("")
    for metrics in ["AUC"] + range(len(EF_FRACTIONS)):
        if metrics == "AUC": print("AUC")
        else: print("EF{}".format(EF_FRACTIONS[metrics]))
        for group in sorted(agg_results):
            for target in agg_results[group]:
                if metrics == "AUC":
                    print("{},{},{:.2f}".format(group, target, np.mean([x["auc"] for x in agg_results[group][target]])))
                else:
                    print("{},{},{:.2f}".format(group, target, np.mean([x["ef"][metrics] for x in agg_results[group][target]])))


def evaluate_logsdirectory(dir):
    """
    The purpose of this function is to merge results of different parameter settings of Bayescreen. It is used
     solely for the purpose of results publication.
    :param dir: Directory to be scanned recursively for logs files  (files ending with .log suffix)
    :return:
    """
    results = []
    results_names = []
    results_file_names = []
    for log in common.find_files_recursively(dir, "*.log"):
        with common.open_file(log) as f:
            results_file_names.append(os.path.basename(log).replace(".log", ""))
            inside = False
            results.append([])
            for line in f:
                if line.startswith("8"): inside = True
                if inside:
                    s_line = line.split(",")
                    results[-1].append(common.to_float(s_line[2]))
                    if len(results_file_names) == 1:
                        results_names.append(s_line[0:2])

    l = len(results[0])
    equal_length = True
    for res in results:
        if len(res) != l:
            equal_length = False
            break
    if not equal_length:
        print ("The result sets differ in size (number of targets")
        for i in range(len(results_file_names)):
            print("{}: {}".format(results_file_names[i], len(results[i])))
        exit()

    print(",," + ",".join(results_file_names))
    for i in range(len(results_names)):
        values = ""
        for res in results:
            values += "{},".format(res[i])
        print("{},{},{}".format(results_names[i][0], results_names[i][1], values.strip(",")))


def main():
    if args.actives is not None and args.inactives is not None:
        print(evaluate_pair(args.actives, args.inactives))
    elif args.directory is not None:
        evaluate_directory(args.directory)
    elif args.logs_directory is not None:
        evaluate_logsdirectory(args.logs_directory)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--actives",
                        metavar='FILE',
                        help="Actives ranking")
    parser.add_argument("-i", "--inactives",
                        metavar='FILE',
                        help="Inactives ranking")
    parser.add_argument("-t", "--type",
                        metavar='[roc|ef]',
                        help="Performance measure -'auc' (area under the ROC curve) "
                             "or 'ef' (enrichment factor)")
    parser.add_argument("-d", "--directory",
                        metavar='DIR',
                        help="Directory with pairs files to evaluate (*actives.out,*inactives.out)")
    parser.add_argument("-l", "--logs-directory",
                        metavar='LOGS_DIR',
                        help="Directory with logs")
    parser.add_argument("--target-index", default="-2", type=int,
                        help="Index of the target name in the directory path")
    parser.add_argument("--group-index", default="-4", type=int,
                        help="Index of the group name in the directory path")
    args = parser.parse_args()

    main()