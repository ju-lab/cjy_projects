import subprocess
import pysam
import os
import sys
import re
from collections import Counter
import numpy as np 


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




reference = '/home/users/team_projects/NormalBRCA/reference/hg19.fa'
ref = pysam.FastaFile(reference)

# Region of interest for capture
#bed = '/home/users/team_projects/NormalBRCA/bed/0813181_Covered.bed'
bed = '/home/users/team_projects/NormalBRCA/analysis/redo.bed'
regions = []
with open(bed, 'r') as f:
    for line in f:
        chromosome, start, end, capture_name = line.strip().split()
        regions.append((Position(chromosome, int(start), int(end)), capture_name))
print(regions)

bams = ['/home/users/team_projects/NormalBRCA/bam/92247-D-29.sorted.whole.postdedup.rg.bam', '/home/users/team_projects/NormalBRCA/bam/114846-D-19.sorted.whole.postdedup.rg.bam']

base_order = dict({'A':0, 'C':1, 'G': 2, 'T':3})

def count_snv_rate(pos, bampath, referencepath):
    basecount = []
    bam = pysam.AlignmentFile(bampath)
    ref = pysam.FastaFile(referencepath)
    
    refbase = ref.fetch(pos.chromosome, pos.start, pos.end).upper()
    counts = bam.count_coverage(pos.chromosome, pos.start, pos.end, quality_threshold=30)

    refCount = counts[base_order[refbase]][0]
    totalCount = counts[0][0] + counts[1][0] + counts[2][0] + counts[3][0]
    return refCount, totalCount, counts[0][0], counts[1][0], counts[2][0], counts[3][0]



def per_position(position, position_name):
    print(position)
    normal_snv_rate_list = []
    tumour_snv_rate_list = []
    with open('vaf_rate_all_capture.txt2', 'a') as f:
        for bam in bams:
            refCount, totalCount, aCount, cCount, gCount, tCount  = count_snv_rate(position, bam, reference)
            f.write(os.path.basename(bam) + '\t' + str(position) + '\t' + str(refCount) + '\t' + str(totalCount) + '\t' + str(aCount) + '\t' + str(cCount) + '\t' + str(gCount) + '\t' + str(tCount) +'\t' + position_name + '\n')

    return 0

if os.path.isfile('vaf_rate_all_capture.txt2'):
    subprocess.call(['rm', '-rf', 'vaf_rate_all_capture.txt2'])

for capture_regions, capture_name in regions:
    for pos in capture_regions:
        print(pos)
        per_position(pos, capture_name)



