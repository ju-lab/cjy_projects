'''Script to identify de novo dominant, autosomal recessive, and comp hets in SSc family
Input file is VCF file that has been left-normalized and decomposed with vt and annotated with VEP. 
2018.06.21 CJY

Different from PDP analysis in that population frequency filter has been loosened to 0.1 from 0.01 
since no variant was found with 0.01 cutoff.

2018.06.27 CJY
Found an important error in determining heterozygous state. Previously only used [0,1] as het.
Now also treats [1,0] as het as well, and also other genotypes. Implemented in function 'variant_type'
'''

import re
import cyvcf2
from collections import Counter

def variant_type(genotype):
    '''convert genotype from cyvcf2 into a string for variant classification'''
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

####################################################################
####################################################################
# COMPOUND HETEROZYGOUS VARIANTS
comphet_candidate_genes = [] # gene name of those that are present as het in affected with at least 2 het variant
for variant in cyvcf2.VCF('/home/users/cjyoon/Projects/rheum/data_processing/01_freebayes/everyone.freebayes.decomposed.norm.vep.centelexcl.vcf.gz'):
    pdp1, pdp2, pdp3, pdp4, ssc1, ssc2, ssc3, ssc4, ssc5 = variant.genotypes
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

    if variant_type(ssc1[0:2]) == 'het' and variant_type(ssc2[0:2]) == 'het' and \
    ((variant_type(ssc3[0:2])=='homo_ref' and variant_type(ssc4[0:2])=='het') or (variant_type(ssc3[0:2])=='het' and variant_type(ssc4[0:2]) == 'homo_ref')) and \
    variant_type(ssc5[0:2])!='homo_alt' and \
    not re.search(r'(intron_variant|regulatory_region_variant|intergenic_variant|downstream|upstream|UTR|non_coding|TF_binding_site_variant|synonymous)', consequence) and (float(max_af_gnomad) < float(freq_threshold)):
        if gene in comphet_candidate_genes:
            # pass
            print(f'{variant.CHROM}\t{variant.POS}\t{variant.REF}\t{variant.ALT[0]}\t{gene}\t{protein_change}\t{sift}\t{polyphen}\t{domains}\t{max_af_gnomad}\t{consequence}') 

        comphet_candidate_genes.append(gene)


# now filter for those variant that are present twice in the list, since comp het needs two different variants
comphet_candidate_genes_count_gt2 = []
for gene, count in Counter(comphet_candidate_genes).items():
    if count >=2:
        comphet_candidate_genes_count_gt2.append(gene)

ncomphet_genes = len(comphet_candidate_genes_count_gt2)

# Now go through the variant list again to write these variants into a file
count = 0 

candidate_file = '/home/users/cjyoon/Projects/rheum/data_processing/01_freebayes/fishing/SSc_candidates_v2.txt'
with open(candidate_file, 'w') as f:        
    f.write('chromosome\tposition\tref\talt\tgene\tprotein_change\tsift\tpolyphen\tdomains\tmax_af_gnomad\tconsequence\tfiltering_strategy\n')
    for variant in cyvcf2.VCF('/home/users/cjyoon/Projects/rheum/data_processing/01_freebayes/everyone.freebayes.decomposed.norm.vep.centelexcl.vcf.gz'):
        pdp1, pdp2, pdp3, pdp4, ssc1, ssc2, ssc3, ssc4, ssc5 = variant.genotypes
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

        if variant_type(ssc1[0:2]) == 'het' and variant_type(ssc2[0:2]) == 'het' and \
        ((variant_type(ssc3[0:2])=='homo_ref' and variant_type(ssc4[0:2])=='het') or (variant_type(ssc3[0:2])=='het' and variant_type(ssc4[0:2]) == 'homo_ref')) and \
        variant_type(ssc5[0:2])!='homo_alt' and not re.search(r'(intron_variant|regulatory_region_variant|intergenic_variant|downstream|upstream|UTR|non_coding|TF_binding_site_variant|synonymous)', consequence) and (float(max_af_gnomad) < float(freq_threshold)):
            if gene in comphet_candidate_genes_count_gt2:
                count += 1
                f.write(f'{variant.CHROM}\t{variant.POS}\t{variant.REF}\t{variant.ALT[0]}\t{gene}\t{protein_change}\t{sift}\t{polyphen}\t{domains}\t{max_af_gnomad}\t{consequence}\tcomp_het\n') 
    
    print(f'variants satisfying compound heterozygous {count} in {ncomphet_genes} genes')

####################################################################
####################################################################
# DE NOVO AUTOSOMAL VARIANT
count = 0
with open(candidate_file, 'a') as f:        
    for variant in cyvcf2.VCF('/home/users/cjyoon/Projects/rheum/data_processing/01_freebayes/everyone.freebayes.decomposed.norm.vep.centelexcl.vcf.gz'):
        pdp1, pdp2, pdp3, pdp4, ssc1, ssc2, ssc3, ssc4, ssc5 = variant.genotypes
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

        if variant_type(ssc1[0:2]) == 'het' and variant_type(ssc2[0:2]) == 'het' and variant_type(ssc3[0:2]) == 'homo_ref' and \
        variant_type(ssc4[0:2]) == 'homo_ref' and variant_type(ssc5[0:2]) == 'homo_ref' and \
        max_af_gnomad < freq_threshold and not re.search(r'(intron_variant|regulatory_region_variant|intergenic_variant|downstream|upstream|UTR|non_coding|TF_binding_site_variant|synonymous)', consequence) :
            f.write(f'{variant.CHROM}\t{variant.POS}\t{variant.REF}\t{variant.ALT[0]}\t{gene}\t{protein_change}\t{sift}\t{polyphen}\t{domains}\t{max_af_gnomad}\t{consequence}\tdenovo_ad\n')
            count += 1

    print(f'variants satisfying autosomal dominant de novo model {count}')


####################################################################
####################################################################
# AUTOSOMAL RECESSIVE VARIANTS
count = 0
with open(candidate_file, 'a') as f:        
    for variant in cyvcf2.VCF('/home/users/cjyoon/Projects/rheum/data_processing/01_freebayes/everyone.freebayes.decomposed.norm.vep.centelexcl.vcf.gz'):
        pdp1, pdp2, pdp3, pdp4, ssc1, ssc2, ssc3, ssc4, ssc5 = variant.genotypes
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

        if variant_type(ssc1[0:2]) == 'homo_alt' and variant_type(ssc2[0:2]) == 'homo_alt' and variant_type(ssc3[0:2]) == 'het' and variant_type(ssc4[0:2]) == 'het' and variant_type(ssc5[0:2]) != 'homo_alt' and max_af_gnomad < freq_threshold and not re.search(r'(intron_variant|regulatory_region_variant|intergenic_variant|downstream|upstream|UTR|non_coding|TF_binding_site_variant|synonymous)', consequence) :
            f.write(f'{variant.CHROM}\t{variant.POS}\t{variant.REF}\t{variant.ALT[0]}\t{gene}\t{protein_change}\t{sift}\t{polyphen}\t{domains}\t{max_af_gnomad}\t{consequence}\tar\n')
            count += 1

    print(f'variants satisfying autosomal recessive {count}')




