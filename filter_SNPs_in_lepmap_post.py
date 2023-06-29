#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
filter_SNPs_in_lepmap_post.py
Version 0.1
Author: Markus Hiltunen
E-mail: markus.hiltunen@ebc.uu.se

Copyright (c) 2020, Johannesson lab
Licensed under the MIT license. See LICENSE file.

Description: filter_SNPs_in_lepmap_post.py is used to filter lepmap posterior
    files using a bed file.
"""

import fileinput

def read_bed(bedfile):
    with open(bedfile, "r") as bed:
        beddict = {}
        for line in bed:
            line = line.strip()
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            chr, start, end = fields[0], int(fields[1]), int(fields[2])
            if chr not in beddict.keys():
                beddict[chr] = []
            beddict[chr].append((start,end))
    return beddict

def main():
    bed = read_bed("/mnt/sda/johannesson_lab/marasmius/marasmius_assemblies/annotations/asm2_repeat_annotation/asm2.2_repeatmasked/asm2.2.fasta.bed")
    for line in fileinput.input():
        printline = True
        line = line.strip()
        if not line.startswith("CHR"):
            fields = line.split("\t")
            chr, pos = fields[0], int(fields[1])
            for region in bed[chr]:
                if pos > region[0] and pos < region[1]:
                    printline = False
                    break
        if printline:
            print(line)


if __name__ == "__main__":
    main()
