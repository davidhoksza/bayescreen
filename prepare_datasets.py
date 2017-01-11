#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Processing fragments' features files in the CSV format. Each line
consists of a fragment SMILES and a list of features. The header lists
features names.
"""

import json
import common
import argparse

__author__ = "David Hoksza"
__email__ = "david.hoksza@mff.cuni.cz"
__license__ = 'X11'


def read_sdf(f):
    molecules = {}
    mol = ""
    mol_id = ""
    for line in f:
        if mol == "":
            mol_id = line.strip("\n")
        mol += line
        if line.startswith("$$$$"):
            molecules[mol_id] = mol
            mol = ""
    return molecules


def main():
    with common.open_file(args.input) as f:
        split_data = json.load(f)
        [fn_inactives, fn_actives] = split_data["files"]

        with common.open_file(args.dir + "/" + fn_actives + ".sdf") as fa:
            molecules = read_sdf(fa)

            fns = [args.dir + "/" + fn_actives + "_train.sdf", args.dir + "/" + fn_actives + "_test.sdf"]
            if args.out_actives_train is not None:
                fns[0] = args.out_actives_train
            if args.out_actives_test is not None:
                fns[1] = args.out_actives_test
            with common.open_file(fns[0], "w") as fta:
                for molecule in split_data["data"]["train"]["ligands"]:
                    fta.write(molecules[molecule["name"]])
            with common.open_file(fns[1], "w") as fta:
                for molecule in split_data["data"]["test"]:
                    if molecule["activity"] == 1:
                        fta.write(molecules[molecule["name"]])

        with common.open_file(args.dir + "/" + fn_inactives + ".sdf") as fi:
            molecules = read_sdf(fi)

            fns = [args.dir + "/" + fn_inactives + "_train.sdf", args.dir + "/" + fn_inactives + "_test.sdf"]
            if args.out_inactives_train is not None:
                fns[0] = args.out_inactives_train
            if args.out_inactives_test is not None:
                fns[1] = args.out_inactives_test
            with common.open_file(fns[0], "w") as fta:
                for molecule in split_data["data"]["train"]["decoys"]:
                    fta.write(molecules[molecule["name"]])
            with common.open_file(fns[1], "w") as fta:
                for molecule in split_data["data"]["test"]:
                    if molecule["activity"] == 0:
                        fta.write(molecules[molecule["name"]])


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input",
                        required=True,
                        help="JSON input file")
    parser.add_argument("-d", "--dir",
                        required=True,
                        help="Directory with sdf files (if different from current dir)",
                        default=".")
    parser.add_argument("-o", "--out-dir",
                        help="Output directory for sdf files (if different from input directory dir)")
    parser.add_argument("--out-actives-train",
                        help="Output file name (if different from the input file names.")
    parser.add_argument("--out-inactives-train",
                        help="Output file name (if different from the input file names.")
    parser.add_argument("--out-actives-test",
                        help="Output file name (if different from the input file names.")
    parser.add_argument("--out-inactives-test",
                        help="Output file name (if different from the input file names.")

    args = parser.parse_args()

    if args.out_dir == "": args.out_dir = args.dir

    main()
