#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description="Convert vcf data to lepmap3 input format.")
parser.add_argument("input_vcf", \
                    help="Input vcf file. Required.", \
                    type = str)
args = parser.parse_args()

def formatGenotype(genotype):
    """Coding genotypes for haploid offspring from a diploid mother,
    following https://sourceforge.net/p/lep-map3/discussion/general/thread/525c83d8fe/?limit=25
    """
    if genotype == "0/0":
        # Coded as AA
        return "1 0 0 0 0 0 0 0 0 0"
    elif genotype == "0/1":
        # Coded as missing data
        return "1"
    elif genotype == "1/1":
        # Coded as AB
        return "0 1 0 0 0 0 0 0 0 0"
    elif genotype == "./.":
        return "1"
    else:
        return "1"


def readVcf(vcffile):
    out = []
    with open(vcffile, "r") as vcf:
        for line in vcf:
            line = line.strip()
            if line.startswith("##"):
                continue
            elif line.startswith("#CHROM"):
                sample_names = "CHR\tPOS\tfemale\tmale\t" + "\t".join(line.split("\t")[9:])
                out.append(sample_names)
            else:
                # Coding parents as mother AB, father AA
                variant_line = "\t".join(line.split("\t")[:2]) + "\t0 1 0 0 0 0 0 0 0 0" + "\t1 0 0 0 0 0 0 0 0 0\t" + "\t".join([ formatGenotype(gen.split(":")[0]) for gen in line.split("\t")[9:]])
                out.append(variant_line)
    return out


def main():
    vcf = readVcf(args.input_vcf)
    with open(args.input_vcf.rstrip(".vcf")+".lepmap.data", "w") as out:
        for l in vcf:
            out.write(l+"\n")

if __name__ == "__main__":
    main()
