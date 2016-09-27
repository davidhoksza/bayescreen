#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Analyzes Bayes activity model.

"""

import argparse
import common
import json
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

__author__ = "David Hoksza"
__email__ = "david.hoksza@mff.cuni.cz"
__license__ = 'X11'

def multipage(filename, figs=None, dpi=200):
    pp = PdfPages(filename)
    if figs is None:
        figs = [plt.figure(n) for n in plt.get_fignums()]
    for fig in figs:
        fig.savefig(pp, format='pdf')
    pp.close()

def analyze(model):
    print("==== General information ====")
    print("Types of fragments: {}".format(model["fragment_types"]))
    print("Tool used for feature values generation: {}".format(model["features_generator"]))
    print("Features: {}".format(",".join(model["features_names"])))

    cnt_bins = model["cnt_bins"]
    probs = model["probabilities"]
    probability_ratios = {}
    for ix_fn in range(len(model["features_names"])):
        fn = model["features_names"][ix_fn]
        probability_ratios[fn] = []
        for ix_bin in range(cnt_bins):
            probability_ratios[fn].append(probs["feature_value_in_actives"][ix_fn][ix_bin] /
                                          probs["feature_value_in_inactives"][ix_fn][ix_bin])

    print("\n==== Relative features importance aggregated over all values of given feature ====")
    sum_ratios = {}
    for fn in probability_ratios:
        sum_ratios[fn] = np.sum(probability_ratios[fn])
    s = np.sum(sum_ratios.values())
    for fn, sr in sorted(sum_ratios.items(), key=lambda x: x[1], reverse=True):
        print("{}:{} ({})".format(fn, s, sr/s))

    mins = model["normalization"]["mins"]
    maxs = model["normalization"]["maxs"]
    print("\n==== Features values importance====")
    print("p(f_i=x|m=A)/p(f_i=x|m=I), feature name, bin")
    all_ratios = []
    for ix_fn in range(len(model["features_names"])):
        fn = model["features_names"][ix_fn]
        bin_size = (maxs[ix_fn] - mins[ix_fn]) / cnt_bins
        for ix_bin in range(cnt_bins):
            all_ratios.append({"value": probability_ratios[fn][ix_bin],
                               "feature": fn,
                               "bin": ix_bin,
                               "range": "({:.3f};{:.3f})".format(mins[ix_fn] + bin_size * ix_bin, mins[ix_fn] + bin_size * (ix_bin + 1))})

    #     fig = plt.figure()
    #     plt.title(fn)
    #     plt.plot(np.linspace(mins[ix_fn], maxs[ix_fn], num=cnt_bins+1), [probability_ratios[fn][0]]+ probability_ratios[fn], drawstyle='steps', fillstyle='full')
    #
    # multipage('multipage1.pdf')

    cnt_printed = 0
    for rec in sorted(all_ratios, key=lambda x: x["value"], reverse=True):
        print("{}, {}, {}, {}".format(rec["value"], rec["feature"], rec["bin"], rec["range"]))
        cnt_printed += 1
        if not args.full and cnt_printed > 20: break;


def main():
    with common.open_file(args.model) as f:
        model = json.load(f)

    if args.log is not None:
        sys.stdout = common.open_file(args.log, "w")
    analyze(model)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-m", "--model",
                        help="Output JSON file containing the resulting model.")
    parser.add_argument("-l", "--log",
                        metavar='FILE',
                        type=str,
                        help="Output log file. If not set, standard output will be used instead.")
    parser.add_argument("-f", "--full",
                        action='store_true',
                        default=False,
                        help="If true, prints all the probabilities. If false, prints only top 20 probabilities")

    args = parser.parse_args()

    main()