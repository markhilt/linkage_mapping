#!/usr/bin/python
# -*- coding: utf-8 -*-

import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("input", help="Input lepmap3 map file", type = str)
parser.add_argument("-n", "--n_markers", help="Only output linkage groups containing -n markers or more [0]", default = 0, type = int)
parser.add_argument("-o", "--output_prefix", help="Output prefix", default = "map", type = str)
args = parser.parse_args()

def readMap(map):
    ''' Read and return linkage map.
    '''
    with open(map, "r") as lmap:
        lgs = {}
        for line in lmap:
            if line.startswith("#***"):
                lg = line.split(" ")[3]
                lgs[lg] = []
            elif not line.startswith("#"):
                fields = line.split("\t")
                tig, coord, gen_pos = fields[0], fields[1], fields[3]
                lgs[lg].append( (tig, coord, gen_pos) )
    return lgs

def main():
    map = readMap(args.input)
    outls = ["Scaffold ID,scaffold position,LG,genetic position"]
    for k, v in map.items():
        if len(v) > args.n_markers:
            for marker in v:
                outls.append("{},{},{},{}".format(marker[0], marker[1], k, marker[2]))

    with open(args.output_prefix+".csv", "w") as out:
        out.write("\n".join(outls))
        out.write("\n")

if __name__ == "__main__":
    main()
