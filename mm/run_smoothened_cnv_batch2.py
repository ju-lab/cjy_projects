import re
import os
import sys
import subprocess
import shlex

def parse_sequenza_output(sequenza_output):
    """extracts aberrant cell fraction and ploidy from sequenza output"""
    with open(sequenza_output, 'r') as f:
        for i, line in enumerate(f):
            if i==2: # read second line (first line is header)
                acf, ploidy, ploidy_mean_cn = line.strip().split()

    return acf, ploidy

def parse_100kb_stat(binning_stat_file):
    """extracts average depth from 100kbcov.covstat file"""
    with open(binning_stat_file, 'r') as f:
        for i, line in enumerate(f):
            if i==1: # read second line (first line is header)
                filename, region, throughput, avg_depth = line.strip().split()
                return avg_depth


def command(tumor_100kb_cov, normal_100kb_cov, tumor_100kb_stat, normal_100kb_stat, sequenza_output):
    """Creates commandline to submit /home/users/sypark/03_Tools/Smoothened_CN/03_report_Tspecific_absolute_CN.py"""    
    acf, ploidy = parse_sequenza_output(sequenza_output)
    tumor_depth = parse_100kb_stat(tumor_100kb_stat)
    normal_depth = parse_100kb_stat(normal_100kb_stat)

    cmd = f"/usr/bin/python /home/users/sypark/03_Tools/Smoothened_CN/03_report_Tspecific_absolute_CN.py {tumor_100kb_cov} {normal_100kb_cov} {tumor_depth} {normal_depth} {acf} {ploidy}"
    print(cmd)
    subprocess.call(shlex.split(cmd))

def main():
    with open('/home/users/cjyoon/Projects/myeloma/sample_info/tumor_match_abs_batch2_20180706.txt', 'r') as f:
        for line in f:
            if not line.startswith('#'):
                sampleID, tumor, normal = line.strip().split()
                tumorID = f'MM01_{sampleID}_T_G'
                normalID = f'MM01_{sampleID}_N_G'
        
                tumor_100kb_cov = f"/home/users/cjyoon/Projects/myeloma/bamfiles/tumor/{tumorID}.bam.mpileup.100kbcov"
                normal_100kb_cov = f"/home/users/cjyoon/Projects/myeloma/bamfiles/normal/{normalID}.bam.mpileup.100kbcov"
                tumor_100kb_stat = f"/home/users/cjyoon/Projects/myeloma/bamfiles/tumor/{tumorID}.bam.mpileup.100kbcov.covstat"
                normal_100kb_stat = f"/home/users/cjyoon/Projects/myeloma/bamfiles/normal/{normalID}.bam.mpileup.100kbcov.covstat"
                sequenza_output = f"/home/users/cjyoon/Projects/myeloma/analysis/sequenza/{sampleID}_confints_CP.txt"

                command(tumor_100kb_cov, normal_100kb_cov, tumor_100kb_stat, normal_100kb_stat, sequenza_output)


if __name__ == '__main__':
    main()
