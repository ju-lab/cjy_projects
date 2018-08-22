'''Script to identify de novo dominant, autosomal recessive, and comp hets in SSc family
Input file is VCF file that has been left-normalized and decomposed with vt and annotated with VEP. 
2018.06.21 CJY

Different from PDP analysis in that population frequency filter has been loosened to 0.1 from 0.01 
since no variant was found with 0.01 cutoff.

Modified 2018.07.01. Fixing variant_type and vaf calculation errors. 
Modiifed 2018.07.02. Output all variants regardless of max_af_vaf and variant classification to include MT variants and just to make the table very inclusive. 
Modified 2018.08.05 v3 Added depth of each sample to allow appropriate depth filtering for better identification of all de novo variants 
'''

import re
import cyvcf2
from collections import Counter
import sys
import numpy as np 
def variant_type(genotype):
    '''convert genotype from cyvcf2 into a string for variant classification
    if variant is not genotyped, then index 0 is -1. need to account for this
    '''
    if genotype[0] == -1:
        return 'NA'
    else:
        if genotype == [0, 0]:
            return 'homo_ref'
        elif genotype == [0, 1] or genotype == [1, 0]:
            return 'het'
        elif genotype != [0, 0] and genotype[0] == genotype[1]:
            return 'homo_alt'
        elif genotype != [0, 1] and genotype != [0, 1] and genotype[0] != genotype[1]:
            return 'het'
        else:
            print('Not prepared for this type of variant genotypes')
            raise ValueError


freq_threshold = 0.1
max_depth_threshold = 100

def find_canonical_annotation(vep_annotation_string): 

    """VEP annotates with many alternative transcripts as well as canonical transcript
    this function finds the canonical transcript within vep_annotation_string.
    If there is no canonical transcript, which is usually the case fore intergenic,
    will just report the first annotation.
    """
    annotations = vep_annotation_string.split(',') 
    return_status = 0
    for annotation in annotations: 
        CANONICAL = annotation.split('|')[26] # CANONICAL
        if CANONICAL == 'YES': 
            return_status = 1
            return annotation 

        if return_status == 0: 
            return vep_annotation_string.split(',')[0]

def calculate_vaf(alt_depth, total_depth):
    '''for a given depth of both total and alt from a variant info, will calulate the vaf'''
    if total_depth != 0:
        return round(float(alt_depth/total_depth), 3)
    else:
        return 'NA'

def sample_vafs(variant):
    pdp1_vaf = calculate_vaf(variant.gt_alt_depths[0], variant.gt_depths[0])
    pdp2_vaf = calculate_vaf(variant.gt_alt_depths[1], variant.gt_depths[1])
    pdp3_vaf = calculate_vaf(variant.gt_alt_depths[2], variant.gt_depths[2])
    pdp4_vaf = calculate_vaf(variant.gt_alt_depths[3], variant.gt_depths[3])
    ssc1_vaf = calculate_vaf(variant.gt_alt_depths[4], variant.gt_depths[4])
    ssc2_vaf = calculate_vaf(variant.gt_alt_depths[5], variant.gt_depths[5])
    ssc3_vaf = calculate_vaf(variant.gt_alt_depths[6], variant.gt_depths[6])
    ssc4_vaf = calculate_vaf(variant.gt_alt_depths[7], variant.gt_depths[7])
    ssc5_vaf = calculate_vaf(variant.gt_alt_depths[8], variant.gt_depths[8])

    return pdp1_vaf, pdp2_vaf, pdp3_vaf, pdp4_vaf, ssc1_vaf, ssc2_vaf, ssc3_vaf, ssc4_vaf, ssc5_vaf

vcf_handle = cyvcf2.VCF('/home/users/cjyoon/Projects/rheum/data_processing/01_freebayes/everyone.freebayes.decomposed.norm.vep.centelexcl.vcf.gz')
writer = cyvcf2.Writer('/home/users/cjyoon/Projects/rheum/data_processing/01_freebayes/everyone.freebayes.decomposed.norm.vep.centelexcl.denovo_v2.vcf', vcf_handle)
for variant in vcf_handle:
    pdp1, pdp2, pdp3, pdp4, ssc1, ssc2, ssc3, ssc4, ssc5 = variant.genotypes
    pdp1_geno = variant_type(pdp1[0:2])
    pdp2_geno = variant_type(pdp2[0:2])
    pdp3_geno = variant_type(pdp3[0:2])
    pdp4_geno = variant_type(pdp4[0:2])
    ssc1_geno = variant_type(ssc1[0:2])
    ssc2_geno = variant_type(ssc2[0:2])
    ssc3_geno = variant_type(ssc3[0:2])
    ssc4_geno = variant_type(ssc4[0:2])
    ssc5_geno = variant_type(ssc5[0:2])

    pdp1_vaf, pdp2_vaf, pdp3_vaf, pdp4_vaf, ssc1_vaf, ssc2_vaf, ssc3_vaf, ssc4_vaf, ssc5_vaf = sample_vafs(variant)
    # print(f'{pdp1_vaf}, {pdp2_vaf}, {pdp3_vaf}, {pdp4_vaf}, {ssc1_vaf}, {ssc2_vaf}, {ssc3_vaf}, {ssc4_vaf}, {ssc5_vaf}')

    pdp1_depth, pdp2_depth, pdp3_depth, pdp4_depth, ssc1_depth, ssc2_depth, ssc3_depth, ssc4_depth, ssc5_depth = variant.gt_depths

    csq = variant.INFO.get('CSQ')
    canonical = find_canonical_annotation(csq)
    consequence = canonical.split('|')[1]
    gene = canonical.split('|')[3]
    protein_change = canonical.split('|')[11]
    sift = canonical.split('|')[36]
    polyphen = canonical.split('|')[37]
    domains = canonical.split('|')[38]

    max_af = find_canonical_annotation(csq).split('|')[57]
    gnomad_af = find_canonical_annotation(csq).split('|')[48]
    if max_af.strip() == '':
        max_af = 0

    if gnomad_af.strip() == '':
        gnomad_af = 0


    max_af_gnomad = max(float(max_af), float(gnomad_af))
    totaldepth = variant.INFO['DP']
    avgdepth = round(float(np.sum(np.array(variant.format('DP'))) / 9), 3)
    # only report rare (10% threshold) that are not intronic/regulatory/synonymous... etc and 
    # not in very high depth region of the genome. 
    #if max_af_gnomad < 0.1 and not re.search(r'(intron_variant|regulatory_region_variant|intergenic_variant|downstream|upstream|UTR|non_coding|TF_binding_site_variant|synonymous)', consequence) and avgdepth < max_depth_threshold:

    if pdp1_vaf != 'NA' and pdp2_vaf != 'NA' and pdp3_vaf != 'NA':
        if len(variant.REF) == 1 and len(variant.ALT[0]) == 1:
            if pdp1_geno == 'het' and pdp2_geno == 'homo_ref' and pdp3_geno == 'homo_ref' and pdp1_vaf > 0.35 and pdp1_vaf < 0.65 and pdp2_vaf < 0.1 and pdp3_vaf < 0.1 and pdp1_depth > 15 and pdp2_depth > 15 and pdp3_depth > 15:
                writer.write_record(variant)

writer.close()
vcf_handle.close()

    


