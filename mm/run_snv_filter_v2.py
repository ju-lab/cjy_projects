'''filter with 8 additional samples
2018.07.10 CJY
'''

import os
import subprocess
import re
import shlex
import pysam
import multiprocessing as mp
filter_directory = '/home/users/cjyoon/Projects/myeloma/analysis/snv_filter_v2'
nin_filter_threshold = 0.01
sample_table = '/home/users/cjyoon/Projects/myeloma/sample_info/mm_sample_table_v2.txt'

def sampleNameBam(bamFile):
    """get @RG SM: information as sample name from BAM header"""
    bam = pysam.AlignmentFile(bamFile)
    name = bam.header['RG'][0]['SM']
    return name


# get list of normal Bams 
normal_bam_list = [] 
with open(sample_table, 'r') as f:
    for line in f:
        if not line.startswith('#'):
            sampleID, tumorBam, normalBam = line.strip().split()
            normal_bam_list.append(normalBam)




def run_filter(sampleID, tumorBam, normalBam):
    tumorID = sampleNameBam(tumorBam)
    normalID = sampleNameBam(normalBam)

    mutect_vcf = f'/home/users/cjyoon/Projects/myeloma/analysis/mutect/{sampleID}.mutect.vcf'
    strelka_vcf = f'/home/users/cjyoon/Projects/myeloma/analysis/strelka/{sampleID}_{tumorID}_{normalID}/results/variants/somatic.snvs.passonly.vcf.gz'
    ok_to_run = 1

    # Check if there's all the input files
    if not os.path.isfile(mutect_vcf):
        print(f"{mutect_vcf} not present")
        ok_to_run = 0

    if not  os.path.isfile(strelka_vcf):
        print(f"{strelka_vcf} not present")
        ok_to_run = 0

    

    # Step 1
    # Apply Not In Normal Filter for Mutect
    cmd = f'python /home/users/cjyoon/scripts/variantFilter/snv/mutect_ninFilter.py -m {mutect_vcf} -s {sampleID} -t {nin_filter_threshold} -o {filter_directory}' 
    print(cmd)
    step1 = subprocess.Popen(shlex.split(cmd))
    step1_status = step1.wait()

    filtered_mutect = os.path.join(filter_directory, re.sub(string=os.path.basename(mutect_vcf), pattern=r'.vcf$', repl='.nin_filter_' + str(nin_filter_threshold) + '.vcf'))
    print(filtered_mutect)

    # check step1 status
    if step1_status == 0:
        pass
    else:
        print(f'{sampleID} failed at step1')
        return -1

    # Step 2
    # Intersect filtered Mutect with Strelka

    cmd = f'python /home/users/cjyoon/scripts/variantFilter/snv/intersect_mutect_strelka.py --mutect {filtered_mutect} --strelka {strelka_vcf} -s {sampleID} -o {filter_directory}'
    print(cmd)
    step2 = subprocess.Popen(shlex.split(cmd))
    step2_status = step2.wait()

    # check step2 status
    if step2_status == 0:
        pass
    else:
        print(f'{sampleID} failed at step2')
        return -1

    intersect_vcf = os.path.join(filter_directory, sampleID + '_mns.vcf') 

    # Step 3a
    # Annotate using Panel of Normal

    normal_bam_list_string = '\t'.join(normal_bam_list)
    cmd = f'python /home/users/cjyoon/scripts/annotate_vcf/pon_annotate.py -i {intersect_vcf} -n {normal_bam_list_string} -o {filter_directory}'
    print(cmd)
    step3a = subprocess.Popen(shlex.split(cmd))
    step3a_status = step3a.wait()

    # check step3a status
    if step3a_status == 0:
        pass
    else:
        print(f'{sampleID} failed at step3a')
        return -1


    # Step 3b
    # Filter using PoN VAF
    pon_vcf = re.sub(r'.vcf$', '.pon.vcf', intersect_vcf)
    cmd = f'python /home/users/cjyoon/scripts/annotate_vcf/pon_filter.py -i {pon_vcf} -t 0.01 -o {filter_directory}'
    print(cmd)
    step3b = subprocess.Popen(shlex.split(cmd))
    step3b_status = step3b.wait()

    # check step3b status
    if step3b_status == 0:
        pass
    else:
        print(f'{sampleID} failed at step3b')
        return -1


    # Step 4
    # Annotate Filtered VCF with VEP
    pon_filtered_vcf = re.sub(r'.vcf$', '.filtered.vcf', pon_vcf)
    cmd = f'python /home/users/cjyoon/scripts/annotate_vcf/annotate_vcf.py -i {pon_filtered_vcf} -o {filter_directory}'
    print(cmd)
    step4 = subprocess.Popen(shlex.split(cmd))
    step4_status = step4.wait()

    # check step3 status
    if step2_status == 0:
        pass
    else:
        print(f'{sampleID} failed at step3')
        return -1

    if step1_status == 0 and step2_status == 0 and step3a_status == 0 and step3b_status == 0 and step4_status == 0:
        return 0
    else:
        return -1

arg_lists = []
with open(sample_table, 'r') as f:
    for line in f:
        if not line.startswith('#'):
            sampleID, tumorBam, normalBam = line.strip().split()
            arg_lists.append((sampleID, tumorBam, normalBam))

# run multicore
with mp.Pool(processes=27) as p:
    for arglist, exitstatus in  zip(arg_lists, p.starmap(run_filter, arg_lists)):
        sample = arglist[0]
        print(f'{sapmle} exit status = {exitstatus}')

