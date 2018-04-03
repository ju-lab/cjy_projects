import subprocess
import shlex
import pysam
import os
import sys
import re
from collections import Counter
import numpy as np 
import argparse
import multiprocessing


class Position():
    ''' python class for handling genomic positions
    0-based
    '''
    def __init__(self, chromosome, start, end, is_bp=None, clipped_reads=None):
        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.is_bp = is_bp
        self.clipped_reads = None

    def __repr__(self):
        return str(self.chromosome + ":" + str(self.start) + '-' + str(self.end))

    def __str__(self):
        return str(self.chromosome + ":" + str(self.start) + '-' + str(self.end))

    def __len__(self):
        return int(self.end - self.start)

    def __iter__(self):
        for i in range(self.start, self.end):
            yield Position(self.chromosome, i, i+1)

    def __next__(self):
        self.start += 1
        self.end += 1
        return self

    def __hash__(self):
        return hash((self.chromosome, self.start, self.end))

    def __eq__(self, other):
        if isinstance(other, Position):
            return self.chromosome == other.chromosome and self.start == other.start and self.end == other.end
        else:
            print("Not of the same class, cannot compare equality")
            return None

    def __lt__(self, other):
        if isinstance(other, Position):
            if self.chromosome < other.chromosome:
                return True
            elif self.chromosome == other.chromosome:
                if self.start < other.start:
                    return True
                else:
                    return False
            else:
                return False

# Order of Bases in Dictionary for easy look up
base_order = dict({'A':0, 'C':1, 'G': 2, 'T':3})

def argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--bam', required=True, help='Input Bam file')
    parser.add_argument('-r', '--reference', required=False, help='Reference Fasta', default='/home/users/team_projects/NormalBRCA/reference/hg19.fa')
    parser.add_argument('-c', '--capture_bed', required=False, help='Capture Bed file', default='/home/users/team_projects/NormalBRCA/bed/0813181_Covered.bed')
    parser.add_argument('--multicore', required=False, default=4, type=int, help='Multicore')


    args = vars(parser.parse_args())
    
    return args['bam'], args['reference'], args['capture_bed'], args['multicore']


def parse_regions_from_bed(bedfile, multicore):
    """read in bedfile and return list of captured region in Position class"""
    regions = dict()
    for i in range(0, multicore):
        regions.update({i:set()})

    
    count = 0 
    with open(bedfile, 'r') as f:
        for line in f:
            chromosome, start, end, capture_name = line.strip().split()
            batchNumber = count % multicore
            regions[batchNumber].add((Position(chromosome, int(start), int(end)), capture_name))
            count+=1 
    return regions


def count_snv_rate(pos, bampath, referencepath):
    basecount = []
    bam = pysam.AlignmentFile(bampath)
    ref = pysam.FastaFile(referencepath)
    
    refbase = ref.fetch(pos.chromosome, pos.start, pos.end).upper()
    counts = bam.count_coverage(pos.chromosome, pos.start, pos.end, quality_threshold=30)

    refCount = counts[base_order[refbase]][0]
    totalCount = counts[0][0] + counts[1][0] + counts[2][0] + counts[3][0]
    return refCount, totalCount, counts[0][0], counts[1][0], counts[2][0], counts[3][0]


def per_position(bam, position, position_name, outputfile, reference):
    print(position)
    normal_snv_rate_list = []
    tumour_snv_rate_list = []
    with open(outputfile, 'a') as f:
        refCount, totalCount, aCount, cCount, gCount, tCount = count_snv_rate(position, bam, reference)
        f.write(os.path.basename(bam) + '\t' + str(position.chromosome) + '\t' + str(position.end) + '\t' + str(refCount) + '\t' + str(totalCount) + '\t' + str(aCount) + '\t' + str(cCount) + '\t' + str(gCount) + '\t' + str(tCount) +'\t' + position_name + '\n')

    return 0


def per_batch(bamfile, region_list, reference, batchNumber):
    """per batch is to utilize multicore"""
    outputfile = os.path.basename(bamfile) + '.count' + str(batchNumber)
    for capture_region, capture_name in region_list:
        for pos in capture_region: # iterate over all positions within list using Position class
            per_position(bamfile, pos, capture_name, outputfile, reference)

    return outputfile

def combine_batchoutputs(bam, outputfiles):
    """combine results from multiple batches into a single file and clean up"""
    outputfile = os.path.basename(bam) + '.count'
    with open(outputfile, 'w') as f:
        for outputbatch in outputfiles:
            with open(outputbatch, 'r') as g:
                for line in g:
                    f.write(line)

            # clean up intermediate batch file
            subprocess.call(shlex.split('rm -rf ' + outputbatch))

    return outputfile

def main():
    # argument parsing
    bam, reference, capture_bed, multicore = argument_parser()

    capture_regions = parse_regions_from_bed(capture_bed, multicore)
    
    # prep argument list for multiprocessing
    arg_list = []
    for i in range(0, multicore):
        arg_list.append([bam, capture_regions[i], reference, i])

    with multiprocessing.Pool(processes=multicore) as pool:
        outputfiles = pool.starmap(per_batch, arg_list)

    print(combine_batchoutputs(bam, outputfiles))



if __name__=='__main__':
    main() 

