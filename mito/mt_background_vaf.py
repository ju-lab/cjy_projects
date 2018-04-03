import pysam
import os
import re
import sys
from collections import Counter
import numpy as np
from multiprocessing import Pool

# input batch number as input
batchNumber = int(sys.argv[1])

# copied chrMT.fa in /Database to here, confirmed identical md5sum to avoid losing this chrMT.fa during sanger maintenance
reference = '/lustre/scratch119/casm/team267ms/ysj/cjy_mt_vaf/chrMT.fa'
ref = pysam.FastaFile(reference)


pcawg_normal_paths = '/lustre/scratch119/casm/team267ms/ysj/cjy_mt_vaf/pcawg_normal_abspath_removed_duplicates.txt'
pcawg_normal_paths_list = []

pcawg_tumour_paths = '/lustre/scratch119/casm/team267ms/ysj/cjy_mt_vaf/pcawg_tumour_abspath_removed_duplicates_removed_size0.txt'
pcawg_tumour_paths_list = []


with open(pcawg_normal_paths, 'r') as f:
    for line in f:
        if not line.startswith('#'): # this removes two duplicated samples 4c0e3a22f60db8def0c849e11239317f.bam 2cad4ccfc32c259f8b4a692191aed827.bam
            normalpaths = line.strip()
            pcawg_normal_paths_list.append(os.path.realpath(normalpaths))


with open(pcawg_tumour_paths, 'r') as f:
    for line in f:
        tumourpaths = line.strip()
        pcawg_tumour_paths_list.append(os.path.realpath(tumourpaths))

print(f'Normals used: {len(pcawg_normal_paths_list)}')
print(f'Tumours used: {len(pcawg_tumour_paths_list)}')

base_order = dict({'A':0, 'C':1, 'G': 2, 'T':3})

def count_snv_rate(mt_pos, bampath, referencepath):
    basecount = []
    bam = pysam.AlignmentFile(bampath)
    ref = pysam.FastaFile(referencepath)
    refbase = ref.fetch('MT', mt_pos-1, mt_pos)
    counts = bam.count_coverage('MT', mt_pos-1, mt_pos, quality_threshold=30)

    refCount = counts[base_order[refbase]][0]
    totalCount = counts[0][0] + counts[1][0] + counts[2][0] + counts[3][0]
    return refCount, totalCount, counts[0][0], counts[1][0], counts[2][0], counts[3][0]



def per_position(position):
    print(position)
    normal_snv_rate_list = []
    tumour_snv_rate_list = []
    with open('./vaf_result3/MT_' + str(position) + '_rate.txt', 'w') as f:
        for normalBam in pcawg_normal_paths_list:
            refCount, totalCount, aCount, cCount, gCount, tCount  = count_snv_rate(position, normalBam, reference)
            f.write('normal\t' + os.path.basename(normalBam) + '\t' + str(position) + '\t' + str(refCount) + '\t' + str(totalCount) + '\t' + str(aCount) + '\t' + str(cCount) + '\t' + str(gCount) + '\t' + str(tCount) + '\n')

        for tumourBam in pcawg_tumour_paths_list:
            refCount, totalCount, aCount, cCount, gCount, tCount = count_snv_rate(position, tumourBam, reference)
            f.write('tumour\t' + os.path.basename(tumourBam) + '\t' + str(position) + '\t' + str(refCount) + '\t' + str(totalCount) + '\t' + str(aCount) + '\t' + str(cCount) + '\t' + str(gCount) + '\t' + str(tCount) + '\n') 

    return 0
    

def batchNumber2mtRegions(batchNumber):
    if batchNumber == 4143:
        return np.arange(16569, 16570)
    elif batchNumber < 4143:
        return np.arange(4*(batchNumber -1)+1, 4*(batchNumber)+1)
    else:
        print('batchNumber cannot be greater than 4143')
        raise ValueError

#mt_positions = batchNumber2mtRegions(batchNumber)

# for rerun just use this line instead of the line above 
mt_positions = [batchNumber]

with Pool(1) as p:
    p.map(per_position, mt_positions)

with open('completed_batch.txt', 'a') as h:
    h.write(str(batchNumber) + '\n')


print('done') 


