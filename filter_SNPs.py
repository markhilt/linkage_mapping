#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
analyze_coding_SNPs.py
Version 0.1
Author: Markus Hiltunen
E-mail: markus.hiltunen@ebc.uu.se

Copyright (c) 2020, Johannesson lab
Licensed under the MIT license. See LICENSE file.

Description: filter_SNPs.py is used to filter vcf files of sites that have reads
    from two variants in haploid samples.
"""

import argparse

parser = argparse.ArgumentParser(description="Filter vcf file.")
parser.add_argument("variants", \
                    help="Variant call file. Required", \
                    type = str)
parser.add_argument("-m", "--mask", \
                    help="Mask heterozygous calls as missing data, instead of \
                    removing them entirely", \
                    action="store_false")
args = parser.parse_args()

def filterAll(fields, f_idx, ad_idx):
    '''Used unless -m.
    '''
    # keep only haploid sites, i.e. sites where all samples have
    # maximum of one read with the minor allele
    printline = True
    for val in fields[f_idx+1:]:
        if val != ".": # Skip missing values
            ad = val.split(":")[ad_idx]
            print(ad)
            if ad != ".":
                ref_n, alt_n = int(ad.split(",")[0]), int(ad.split(",")[1])
                if ref_n > 1 and alt_n > 1:
                    printline = False
                    break
    if printline:
        print("\t".join(fields))

def mask(fields, f_idx, ad_idx):
    ''' If user chooses to mask heterozygous samples instead of filtering
    those lines completely.
    '''
    parsedline = fields
    for idx, val in enumerate(fields):
        if idx > f_idx and val != ".":
            ad = val.split(":")[ad_idx]
            if ad != ".":
                ref_n, alt_n = int(ad.split(",")[0]), int(ad.split(",")[1])
                if ref_n > 1 and alt_n > 1:
                    parsedline[idx] = "."
    print("\t".join(parsedline))

def main():
    with open(args.variants, "r") as vcf:
        for line in vcf:
            line = line.strip()
            if line.startswith("#"):
                print(line)
                if line.startswith("#CHROM"):
                    # FInd in which field samples start
                    fields = line.split("\t")
                    for f_idx, f in enumerate(fields):
                        if f == "FORMAT":
                            break

            else:
                fields = line.split("\t")
                #sample_values = fields[f_idx+1:]
                printline = True

                # find which field is allelic depth
                for ad_idx, tag in enumerate(fields[f_idx].split(":")):
                    if tag == "AD":
                        break
                if args.mask:
                    filterAll(fields, f_idx, ad_idx)

                else:
                    mask(fields, f_idx, ad_idx)

if __name__ == "__main__":
    main()
