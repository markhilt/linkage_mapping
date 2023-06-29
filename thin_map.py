#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
thin_map.py
Version 0.1
Author: Markus Hiltunen
E-mail: markus.hiltunen@ebc.uu.se

Copyright (c) 2020, Johannesson lab
Licensed under the MIT license. See LICENSE file.

Description: thin_map.py is used to collapse markers that have the same cM
    distance in a genetic map. Assumes the position is in column 2 of a tab
    delimited text file.
"""

import argparse

parser = argparse.ArgumentParser(description="Filter genetic map.")
parser.add_argument("map", \
                    help="Map file. Required", \
                    type = str)
args = parser.parse_args()

def main():
    with open(args.map, "r") as map:
        printed_positions = []
        for line in map:
            line = line.strip()
            if line.startswith("#"):
                print(line)
            else:
                fields = line.split("\t")
                recomb_pos = fields[0] + fields[1]
                if recomb_pos not in printed_positions:
                    print(line)
                    printed_positions.append(recomb_pos)

if __name__ == "__main__":
    main()
