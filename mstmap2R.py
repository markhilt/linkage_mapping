#!/usr/bin/python
# -*- coding: utf-8 -*-

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("mstmap_output", help="Input mstmap output file", type = str)
parser.add_argument("-n", "--n_markers", help="Only output linkage groups containing -n markers or more [0]", default = 0, type = int)
args = parser.parse_args()

def main():
    with open(args.mstmap_output, "r") as mstmap:
        linkage_groups = {}
        markers = []
        n_markers = 0
        counter = 0
        insidegroup = False
        for line in mstmap:
            line = line.strip()
            if line.startswith("group"):
                insidegroup = True
                groupname = line.split(" ")[1]
                linkage_groups[groupname] = []
                continue

            elif line.startswith(";ENDOFGROUP"):
                insidegroup = False

            if insidegroup == True and line[0] != ";":
                locus = line.split("\t")[0]
                position = line.split("\t")[1]
                print(groupname, position, locus)

                linkage_groups[groupname].append( (position, locus) )

        outfile = args.mstmap_output.split("/")[-1].rstrip(".txt")+".table.txt"
        with open(outfile, "w") as out:
            out.write("group\tposition\tlocus\n")
            for key, value in linkage_groups.items():
                if len(value) > args.n_markers:
                    for marker in value:
                        out.write("{}\t{}\t{}\n".format(key, marker[0], marker[1]))

if __name__ == "__main__":
    main()
