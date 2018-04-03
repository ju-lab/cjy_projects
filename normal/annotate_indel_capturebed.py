"""
This code adds extra column to the count matrix
1. whether a position is within 10bp of a known indel region called from Varscan2 indel caller
2. whether a position is within the capture bed
April 3 2018 Jongsoo Yoon (cjyoon@kaist.ac.kr)
"""
import cyvcf2
import pandas as pd
import numpy as np 
import re 
import os

# import count data
inputfile = '../analysis/normal_brca_q30_removeclip.count_with_af.txt'
countData = pd.read_table(inputfile, sep='\t')



def indel_positions(indel_vcf):
    indel_positions = dict()
    for chrom in ['chr' + str(i) for i in range(1,23)] + ['chrX', 'chrY', 'chrM']:
        indel_positions.update({chrom:set()})
        
    for variant in cyvcf2.VCF(indel_vcf):
        indel_positions[variant.CHROM].add(variant.POS)
    
    return indel_positions


def get_capture_regions(capture_bed):
    capture_positions = dict()
    for chrom in ['chr' + str(i) for i in range(1,23)] + ['chrX', 'chrY', 'chrM']:
        capture_positions.update({chrom:set()}) 
        
    with open(capture_bed, 'r') as f:
        for line in f:
            chromosome, start, end, *args = line.strip().split()
            capture_positions[chromosome].add((int(start), int(end)))
            
    return capture_positions


# import capture bed file
capture_bed = '/home/users/team_projects/NormalBRCA/bed/0813181_Covered.bed'
capture_regions = get_capture_regions(capture_bed)

capture_df = pd.read_table(capture_bed, header=None, delim_whitespace=True)

# import indel for each bams
indel_1 = indel_positions('../variant_call/varscan_indel/114846-D-19.sorted.whole.postdedup.rg.q30.removeclip.indel.vcf.gz')
indel_9 = indel_positions('../variant_call/varscan_indel/92247-D-29.sorted.whole.postdedup.rg.q30.removeclip.indel.vcf.gz')


def is_in(row, pos):
    """vectorized function to check whether a given position is between the 'start' column and 'end' column"""
    return (row['start'] - pos) * (row['end'] - pos) <=0

def is_in_capture(row):
    """vectorized function to check whether a given position (specified by 'chromosome' and 'position' column)
    are in a given capture bed
    """
    df_bed = capture_df.iloc[:, 0:3]
    df_bed.columns = ['chromosome', 'start', 'end']
    
    df_chrom = df_bed.loc[df_bed['chromosome'] == row['chromosome']]
    truth_value = sum(df_chrom.apply(is_in, pos=row['position'], axis=1))
    
    if truth_value >0:
        return True
    else:
        return False
    



def is_close_to_indel(row, distance_threshold = 10):
    """Vectorized function to check whether a position (specified by 'chromosome' and 'position' column)
    is close to an indel. 
    Return value will be True if it is close to indel locations
    """
    if row['bam'] == '114846-D-19.sorted.whole.postdedup.rg.q30.removeclip.bam':
        indels = indel_1
    elif row['bam'] == '92247-D-29.sorted.whole.postdedup.rg.q30.removeclip.bam':
        indels = indel_9
    
    # now check the distance to indels within the same chromosome
    close_to_indel_count = np.sum(abs(np.array(list(indels[row['chromosome']]))- (row['position'])) < distance_threshold  )
    
    if close_to_indel_count > 0 :
        return True
    else:
        return False


countData['close_to_indel'] = countData.apply(is_close_to_indel, axis=1)
countData['is_in_capture'] = countData.apply(is_in_capture, axis=1)

outputpath = re.sub(string=inputfile, pattern=r'.txt$', repl='.indel_bed_annotated.txt')

countData.to_csv(outputpath, sep='\t', index=False)
