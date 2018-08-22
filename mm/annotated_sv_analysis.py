'''This script reads in all annotated SV files and tries to identify common SVs
2018.07.16 CJY
'''
from collections import Counter

samples = [] 
with open('/home/users/cjyoon/Projects/myeloma/sample_info/mm_sample_table_v2.txt') as f:
	for line in f:
		if not line.startswith('#'):
			sampleID, tumorBam, normalBam = line.strip().split()
			samples.append(sampleID)

fusion_gene_pairs = []

def annotation_split(bp_anno):
	geneName, strand, exon, intron, variant_class, nearest = bp_anno.split('|')
	return geneName, strand, exon, intron, variant_class, nearest

def extract_fusion_genes(bp1_anno, bp2_anno):
	bp1_geneName, bp1_strand, bp1_exon, bp1_intron, bp1_variant_class, bp1_nearest = annotation_split(bp1_anno)
	bp2_geneName, bp2_strand, bp2_exon, bp2_intron, bp2_variant_class, bp2_nearest = annotation_split(bp2_anno)

	if bp1_geneName == '':
		bp1_geneName = bp1_nearest

	if bp2_geneName == '':
		bp2_geneName = bp2_nearest

	fusion_pair = [bp1_geneName, bp2_geneName]
	fusion_pair.sort()
	return fusion_pair

def extract_exonic_fusion_genes(bp1_anno, bp2_anno):
	bp1_geneName, bp1_strand, bp1_exon, bp1_intron, bp1_variant_class, bp1_nearest = annotation_split(bp1_anno)
	bp2_geneName, bp2_strand, bp2_exon, bp2_intron, bp2_variant_class, bp2_nearest = annotation_split(bp2_anno)
	if bp1_geneName == '' or bp2_geneName == '':
		return ['NA', 'NA']
	elif bp1_variant_class in ['intron_variant', 'exonic_variant'] and bp2_variant_class in ['intron_variant', 'exonic_variant']:
		fusion_pair = [bp1_geneName, bp2_geneName]
		fusion_pair.sort()
		return fusion_pair
	else:
		return ['NA', 'NA']


fusion_gene_pairs_all = []
for sampleID in samples:
	fusion_gene_pair = set()
	annotated_sv = f'/home/users/cjyoon/Projects/myeloma/analysis/sv_intersect/{sampleID}_sv_filtered.intersectFilter.vcf.svanno.txt'
	with open(annotated_sv, 'r') as f:
		for line in f:
			bp1, bp2, svtype, orientation, bp1_anno, bp2_anno = line.strip().split()
			a, b = (extract_fusion_genes(bp1_anno, bp2_anno))
			fusion_gene_pair.add((a, b))
	fusion_gene_pairs_all += list(fusion_gene_pair)

# print(Counter(fusion_gene_pairs_all))


fusion_exonic_gene_pairs_all = []
for sampleID in samples:
	fusion_exonic_gene_pair = set()
	annotated_sv = f'/home/users/cjyoon/Projects/myeloma/analysis/sv_intersect/{sampleID}_sv_filtered.intersectFilter.vcf.svanno.txt'
	with open(annotated_sv, 'r') as f:
		for line in f:
			bp1, bp2, svtype, orientation, bp1_anno, bp2_anno = line.strip().split()
			a, b = (extract_exonic_fusion_genes(bp1_anno, bp2_anno))
			if a!= 'NA' and b!='NA':
				fusion_exonic_gene_pair.add((a, b))
	fusion_exonic_gene_pairs_all += list(fusion_exonic_gene_pair)

print(Counter(fusion_exonic_gene_pairs_all))
