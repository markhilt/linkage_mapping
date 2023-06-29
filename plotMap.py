#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
plotMap.py
Version 0.1
Author: Markus Hiltunen
E-mail: markus.hiltunen@ebc.uu.se

Copyright (c) 2020, Johannesson lab
Licensed under the MIT license. See LICENSE file.

Description: plotMap.py is used to plot a linkage map to SVG format.
"""

import argparse
import math
import drawSvg as draw

parser = argparse.ArgumentParser(description="Plot linkage map.")
parser.add_argument("map", \
                    help="Map file. Tab delimited with the header: group position locus. \
                    Required", \
                    type = str)
parser.add_argument("fai", \
                    help="Fasta index file. Required for physical distances", \
                    type = str)
parser.add_argument("-r", "--reverse_LGs", help="Comma-separated list of linkage \
                    groups to reverse the marker order of", type = str)
args = parser.parse_args()


class Marker:
    def __init__(self, chr_name, phys_pos, gen_pos):
        self.chr_name = str(chr_name) # Name of the chromosome
        self.phys_pos = int(phys_pos) # Physical coordinate
        self.gen_pos = float(gen_pos) # Genetic position

    def __str__(self):
        return "Marker: {0}\n\t\
        Physical coordinate: {1} bp\n\t\
        Genetic position: {2}".format(  self.chr_name, \
                                        str(self.phys_pos),
                                        str(self.gen_pos))

class Chromosome:
    def __init__(self, name):
        self.lg_name = name # Name of the linkage group
        self.chr_name = "" # Name of the corresponding chromosome
        self.markers = []  # List of marker objects
        self.phys_start = 0 # Physical starting coordinate, always 0
        self.phys_end = 0 # Physical ending coordinate
        self.gen_start = 0 # Genetic starting coordinate in cM
        self.gen_end = 0 # Genetic ending coordinate in cM

    def __str__(self):
        return "Chromosome: {0}\n\t\
            LG: {1}\n\t\
            Length: {2} bp\n\t\
            N markers: {3}\n\t\
            Total distance: {4} cM".format( self.chr_name, \
                                                self.lg_name, \
                                            str(self.phys_end), \
                                            len(self.markers), \
                                            str(self.gen_end))
    def addMarker(self, chr_name, phys_pos, gen_pos):
        self.markers.append(Marker(chr_name, phys_pos, gen_pos)) # Add marker

    def nameChromosome(self):
        '''Find the most common chromosome name in the markers.
        '''
        from collections import Counter
        marker_chr_names = Counter( [marker.chr_name for marker in self.markers] )
        self.chr_name = marker_chr_names.most_common(1)[0][0]

    def updateCoordinates(self, physical_length):
        '''Update physical and genetic positions of the chromosome. For physical
        position, use the given number that's extracted from a fasta index file.
        For genetic position, find the highest and lowest numbers that are on the
        correct chromosome.
        '''
        if physical_length:
            self.phys_end = int(physical_length)
        else:
            self.phys_end = max([marker.phys_pos for marker in self.markers \
                                if marker.chr_name == self.chr_name])
        gen_pos_lst = [ marker.gen_pos for marker in self.markers \
                        if marker.chr_name == self.chr_name]
        self.gen_start = min(gen_pos_lst)
        self.gen_end = max(gen_pos_lst)

    def reverseMarkers(self):
        '''Reverse the order of the markers
        '''
        self.markers = [Marker( marker.chr_name, \
                                marker.phys_pos, \
                                self.gen_end - marker.gen_pos) \
                        for marker in self.markers]

def readFai(fai_file):
    '''Generate fasta lengths
    '''
    with open(fai_file, "r") as fai:
        for line in fai:
            line = line.strip()
            fields = line.split("\t")
            yield fields[0], int(fields[1])

def draw_scalebar(sb_len, scaling_factor, x, y, unit, tick_dist, side):
    '''Return a drawable scalebar.
    '''
    scalebar = [] # Initiate the scale bar
    scalebar.append(draw.Line(x, y, x, y-sb_len, stroke='black', stroke_width=1.5, fill='none'))

    # Add the unit at the center of the bar, at the side given.
    unit_y_pos = y-sb_len/2
    unit_x_pos = x + 40 if side == "right" else x - 40
    skew = 'skewY(-90)' if side == "right" else 'skewY(90)'
    scalebar.append(draw.Text(unit, 15, unit_x_pos, unit_y_pos, center=0.6, \
                    text_anchor="middle", fill='black'))#, transform = skew)

    # Add tick marks.
    tick_y_pos = y
    x_target = x+5 if side == "right" else x-5
    n = 0
    tick_text_x_pos = x + 15 if side == "right" else x - 15
    while tick_y_pos > y-sb_len:
        scalebar.append(draw.Line(x, tick_y_pos, x_target, tick_y_pos, \
                        stroke='black', stroke_width=1, fill='none'))
        # Draw the number at the tick mark
        scalebar.append(draw.Text(str(n), 10, tick_text_x_pos, tick_y_pos, \
                        center=0.6, text_anchor="middle", fill='black'))
        # We want tick marks every tick_dist untits
        tick_y_pos -= tick_dist * scaling_factor
        n += tick_dist
    return scalebar

def draw_chromosome(chr, lg_sf, chrom_sf, x, y):
    '''Return a drawable LG and chromosome.
    '''
    chromosome = []

    # Draw LG as vertical line.
    scaled_gen_length = chr.gen_end * lg_sf
    chromosome.append(draw.Line(x, y, x, y-scaled_gen_length, stroke='black', \
                                stroke_width=1, fill='none'))

    # Draw chromosome. Two parallel lines, connected at top and bottom with
    # arcs.
    chrom_x_pos = x + 30
    scaled_phys_length = (chr.phys_end / 1000000) * chrom_sf # First recount to Mb
    chrom_width = 8
    chromosome.append(  draw.Line(chrom_x_pos-chrom_width/2, y, \
                        chrom_x_pos-chrom_width/2, y-scaled_phys_length, \
                        stroke='black', stroke_width=0.7, fill='none'))
    chromosome.append(  draw.Line(chrom_x_pos+chrom_width/2, y, \
                        chrom_x_pos+chrom_width/2, y-scaled_phys_length, \
                        stroke='black', stroke_width=0.7, fill='none'))
    # Arcs to connect the two lines
    chromosome.append(  draw.Arc(chrom_x_pos, y, \
                        chrom_width/2, 0, 180, cw=False, \
                        stroke='black', stroke_width=0.7, fill='none'))
    chromosome.append(  draw.Arc(chrom_x_pos, y-scaled_phys_length, \
                        chrom_width/2, 180, 360, cw=False, \
                        stroke='black', stroke_width=0.7, fill='none'))

    # Add tick marks on the LG for the marker positions, and connect LG positions
    # chromosome positions. Draw a red tick mark if the chrom_pos is on another
    # chromosome, otherwise black. If order is correct, draw the connection in
    # green, otherwise red.
    for marker in chr.markers:
        # Calculate LG and chromosome positions
        lg_pos = y - marker.gen_pos * lg_sf
        chr_pos = y - (marker.phys_pos/1000000) * chrom_sf
        # Draw LG position
        col = "black" if chr.chr_name == marker.chr_name else "magenta"
        chromosome.append(draw.Line(x-2.5, lg_pos, \
                                    x+2.5, lg_pos, stroke=col, \
                                    stroke_width=0.5, fill='none'))
        # Draw chromosome position, if we're on the right chromosome. ALso
        # draw a connecting line between the LG and chromosome.
        # Else don't draw anything here.
        if col == "black":
            chromosome.append(draw.Line(chrom_x_pos-chrom_width/2, chr_pos, \
                                        chrom_x_pos+chrom_width/2, chr_pos, \
                                        stroke=col, \
                                        stroke_width=0.5, fill='none'))
            chromosome.append(draw.Line(x+2.5, lg_pos, \
                                        chrom_x_pos-chrom_width/2, chr_pos, \
                                        stroke="green", \
                                        stroke_width=0.1, fill='none'))


    # Draw a heatmap of marker density

    # Add titles for the LG and chromosome
    label_x_pos = x + 15
    label_y_pos = y + 25
    chromosome.append(  draw.Text(chr.lg_name, 20, label_x_pos, label_y_pos, center=0.6, \
                        text_anchor="middle", fill='black'))

    return chromosome

def draw_linkage_map(map):
    '''Draw the linkage map.
    '''
    # First define drawing coordinate system. Will depend on number and lengths
    # of linkage groups
    x = 100 + 100 * len(map.keys()) # 100 times the number of LGs + 100 for scale bars and borders
    y = x / 2 # I want a roughly 1:2 scaling between x and y, to fit across the top half of a journal page.

    # Find scaling factors for cM and Mb based on the longest LG and sequence
    cM_max = max( [chr.gen_end for chr in map.values()] )
    Mb_max = max( [chr.phys_end for chr in map.values()] ) / 1000000
    # Leave 20% y margins: 15% on top, 5% on bottom
    cM_scaling_factor = (0.8 * y) / cM_max
    Mb_scaling_factor = (0.8 * y) / Mb_max

    # Initiate drawing board
    d = draw.Drawing(x, y, origin = (0,0), displayInline=False)
    x_pos, y_pos = 55, 0.90*y # Define cursor starting position

    # Draw a cM scale bar on the left side. 10 cM between tick marks
    scalebar_max = math.ceil(cM_max * cM_scaling_factor)
    d.extend(draw_scalebar( scalebar_max, cM_scaling_factor, \
                            x_pos, y_pos, "cM", 10, "left"))

    # Draw chromosomes. First sort them.
    x_pos += 30
    chrom_order = sorted(list(map.keys()))
    for chrom in chrom_order:
        print("Drawing:")
        print(map[chrom])
        d.extend(draw_chromosome(   map[chrom], cM_scaling_factor, \
                                    Mb_scaling_factor, x_pos, y_pos))
        x_pos += 90

    # Draw a Mb scale bar on the right side. 0.5 Mb between tick marks
    scalebar_max = math.ceil(cM_max * cM_scaling_factor)
    d.extend(draw_scalebar( scalebar_max, Mb_scaling_factor, \
                            x_pos, y_pos, "Mb", 0.5, "right"))

    d.setPixelScale(2)
    d.saveSvg('linkage-map.svg')

def readMap(mapfile, faifile, revstr):
    '''Read the map and store in memory as a dict of chromosomes.
    '''
    map_dict = {}
    with open(mapfile, "r") as map:
        print("Reading " + mapfile)
        for line in map:
            line = line.strip()
            if line.startswith("group"):
                continue
            fields = line.split("\t")
            LG_name, gen_pos, marker = fields[0], float(fields[1]), fields[2]
            chr_name, phys_pos = marker.split("_")[0], int(marker.split("_")[-1])
            if LG_name not in map_dict.keys():
                map_dict[LG_name] = Chromosome(LG_name)
            map_dict[LG_name].addMarker(chr_name, phys_pos, gen_pos)

    print("Updating coordinates")
    fasta_lengths = dict(readFai(faifile))
    reverse_these = revstr.split(",")
    for chr_name in map_dict.keys():
        map_dict[chr_name].nameChromosome()
    for chrom in map_dict.values():
        chrom.updateCoordinates(fasta_lengths[chrom.chr_name])
        if chrom.lg_name in reverse_these:
            map_dict[chrom.lg_name].reverseMarkers()
    return map_dict

def printStats(map):
    '''Calculate some statistics from the linkage map.
    '''
    # Map distance stats
    gen_lengths = [chr.gen_end for chr in map.values()]
    max_cM, min_cM = max(gen_lengths), min(gen_lengths)
    tot_cM = sum(gen_lengths)
    avg_cM = tot_cM / len(gen_lengths)
    # Marker stats
    marker_numbers = [len(chr.markers) for chr in map.values()]
    max_markers, min_markers = max(marker_numbers), min(marker_numbers)
    tot_markers = sum(marker_numbers)
    avg_markers = tot_markers / len(marker_numbers)
    print(  "Map statistics:\n\t\
            Total distance: {} cM\n\t\
            Maximum LG distance: {}\n\t\
            Minimum LG distance: {}\n\t\
            Average distance per LG: {}\n\t\
            Total # markers: {}\n\t\
            Maximum # markers for one LG: {}\n\t\
            Minimum # markers for one LG: {}\n\t\
            Average # markers per LG: {}\n\t\
            ".format(   tot_cM, max_cM, min_cM, avg_cM, \
                        tot_markers, max_markers, min_markers, avg_markers))


def main():
    print("Reading map")
    map_dict = readMap(args.map, args.fai, args.reverse_LGs)
    print("Drawing map")
    draw_linkage_map(map_dict)
    printStats(map_dict)

if __name__ == "__main__":
    main()
