#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("corr", help="Input corr file", type = str)
parser.add_argument("faidx", help="Input faidx file", type = str)
#parser.add_argument("-p", "--p_value", help="p-value cutoff for chi2 test [0.05]", default = 0.05, type = float)
args = parser.parse_args()

def readFaidx(input):
    faidx = {}
    with open(input, "r") as f:
        for line in f:
            line = line.strip()
            tig = line.split("\t")[0]
            length = int(line.split("\t")[1])
            faidx[tig] = length
    return faidx

def readCorr(input):
    corr = {}
    with open(input, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                tig_name = line.split(" ")[0][1:]
            else:
                tig_order = line.split(" ")
                corr[tig_name] = tig_order
    return corr

def calcPositions(faidx, corr):
    bedpos = {}
    for k,v in corr.items():
        bedpos[k] = []
        pos = 1 # Bed coordinates start at 1
        for merged_tig in v:
            merged_tig_name = merged_tig[:-1]
            merged_tig_ori = merged_tig[-1]
            tig_len = faidx[merged_tig_name]
            pos += tig_len
            if merged_tig_ori == "?":
                bedpos[k].append((pos-tig_len, pos))
            pos += 100 # Gap size is 100
    return bedpos

def writeBed(bedpositions):
    with open("converted.bed", "w") as out:
        for k,v in bedpositions.items():
            for i in v:
                out.write("{}\t{}\t{}\n".format(k, i[0], i[1]))

def main():
    faidx = readFaidx(args.faidx)
    corr = readCorr(args.corr)
    bedpositions = calcPositions(faidx, corr)
    writeBed(bedpositions)


if __name__ == "__main__":
    main()
