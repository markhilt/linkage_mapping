#!/usr/bin/python
# -*- coding: utf-8 -*-

import argparse
from scipy import stats

parser = argparse.ArgumentParser()
parser.add_argument("vcf", help="Input vcf file", type = str)
parser.add_argument("-p", "--p_value", help="p-value cutoff for chi2 test [0.05]", default = 0.05, type = float)
args = parser.parse_args()


def readVcf(vcf):
    vcf_out = {}
    with open(vcf, "r") as vcf_file:
        for line in vcf_file:
            line = line.strip()
            if line.startswith("#"):
                continue
            else:
                fields = line.split("\t")
                marker = fields[0] + "_" + fields[1]
                genotypes = [int(gt.split(":")[0]) for gt in fields[9:] if gt[0] != "."]
                ch2 = stats.chisquare(genotypes)
                if ch2.pvalue > args.p_value:
                    print(line)

                n_ref = genotypes.count(0)
                n_alt = genotypes.count(1)

                #print(n_ref,n_alt)

def main():
    readVcf(args.vcf)

if __name__ == "__main__":
    main()
