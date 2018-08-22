# Summarize SNVs in Multiple Myeloma WGS (n=29)
# Jongsoo Yoon (cjyoon@kaist.ac.kr)
# Modified from `summarize_filtered_snv.R` to take output from `~/scripts/annotate_vcf/tidy_annotated_vcf.py` that also extracts
# VAF info from 'mutect' based VCF format.
# April 10 2018
# 2018.07.15 modified from v2. Run with additional 8 samples (total 37) after PON filters

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggbio)
library(cowplot)
library(stringr)

df_4638 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/tidy_data_v2/4638_mns.pon.filtered.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings=''))
df_5685 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/tidy_data_v2/5685_mns.pon.filtered.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings=''))
df_3897 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/tidy_data_v2/3897_mns.pon.filtered.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings=''))
df_4663 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/tidy_data_v2/4663_mns.pon.filtered.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings=''))
df_4031 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/tidy_data_v2/4031_mns.pon.filtered.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings=''))
df_4132 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/tidy_data_v2/4132_mns.pon.filtered.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings=''))
df_4043 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/tidy_data_v2/4043_mns.pon.filtered.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings=''))
df_3869 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/tidy_data_v2/3869_mns.pon.filtered.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings=''))
df_4491 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/tidy_data_v2/4491_mns.pon.filtered.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings=''))
df_3346 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/tidy_data_v2/3346_mns.pon.filtered.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings=''))
df_4509 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/tidy_data_v2/4509_mns.pon.filtered.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings=''))
df_5702 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/tidy_data_v2/5702_mns.pon.filtered.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings=''))
df_4644 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/tidy_data_v2/4644_mns.pon.filtered.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings=''))
df_4527 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/tidy_data_v2/4527_mns.pon.filtered.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings=''))
df_4489 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/tidy_data_v2/4489_mns.pon.filtered.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings=''))
df_4515 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/tidy_data_v2/4515_mns.pon.filtered.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings=''))
df_2643 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/tidy_data_v2/2643_mns.pon.filtered.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings=''))
df_4593 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/tidy_data_v2/4593_mns.pon.filtered.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings=''))
df_3796 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/tidy_data_v2/3796_mns.pon.filtered.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings=''))
df_2581 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/tidy_data_v2/2581_mns.pon.filtered.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings=''))
df_5635 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/tidy_data_v2/5635_mns.pon.filtered.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings=''))
df_2665 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/tidy_data_v2/2665_mns.pon.filtered.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings=''))
df_2972 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/tidy_data_v2/2972_mns.pon.filtered.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings=''))
df_3396 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/tidy_data_v2/3396_mns.pon.filtered.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings=''))
df_3018 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/tidy_data_v2/3018_mns.pon.filtered.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings=''))
df_2616 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/tidy_data_v2/2616_mns.pon.filtered.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings=''))
df_3971 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/tidy_data_v2/3971_mns.pon.filtered.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings=''))
df_5032 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/tidy_data_v2/5032_mns.pon.filtered.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings=''))
df_3083 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/tidy_data_v2/3083_mns.pon.filtered.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings=''))
df_2685 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/tidy_data_v2/2685_mns.pon.filtered.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings=''))
df_3373 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/tidy_data_v2/3373_mns.pon.filtered.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings=''))
df_2655 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/tidy_data_v2/2655_mns.pon.filtered.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings=''))
df_4730 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/tidy_data_v2/4730_mns.pon.filtered.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings=''))
df_5580 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/tidy_data_v2/5580_mns.pon.filtered.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings=''))
df_4529 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/tidy_data_v2/4529_mns.pon.filtered.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings=''))
df_5650 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/tidy_data_v2/5650_mns.pon.filtered.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings=''))
df_2709 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/tidy_data_v2/2709_mns.pon.filtered.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings=''))
df_all = rbind(df_4638,df_5685,df_3897,df_4663,df_4031,df_4132,df_4043,df_3869,df_4491,df_3346,df_4509,df_5702,df_4644,df_4527,df_4489,df_4515,df_2643,df_4593,df_3796,df_2581,df_5635,df_2665,df_2972,df_3396,df_3018,df_2616,df_3971,df_5032,df_3083,df_2685,df_3373,df_2655,df_4730,df_5580,df_4529,df_5650,df_2709)
df_all = df_all %>% mutate(sampleName = as.character(sampleName))
View(df_all)
########################################
# VAF of all mutations called in a histogram 
########################################
ggplot(df_all, aes(VAF)) + geom_histogram()


# Landscape plot

########################################
# Number of SNVs in Exons and WGS
########################################

noncoding_mutation_count = df_all %>% select(SYMBOL, sampleName) %>% filter(is.na(SYMBOL)) %>% group_by(sampleName) %>% tally()
noncoding_mutation_count$type='Non-coding'
coding_mutation_count = df_all %>% select(SYMBOL, sampleName) %>% filter(!(is.na(SYMBOL))) %>% group_by(sampleName) %>% tally() 
coding_mutation_count$type='Coding'

combined_count = rbind(noncoding_mutation_count, coding_mutation_count)

# Median value 
snv_count_bysample = df_all %>% group_by(sampleName) %>% tally() 
median(snv_count_bysample$n) # median snv =  5818

df_all %>% group_by(sampleName) %>% tally() 

########################################
# Import and summarize indel data
########################################
# import pysam
# def sampleNameBam(bamFile):
#   """get @RG SM: information as sample name from BAM header"""
# bam = pysam.AlignmentFile(bamFile)
# name = bam.header['RG'][0]['SM']
# return name
# 
# sample_list = []
# with open('mm_sample_table_v2.txt', 'r') as f:
#   for line in f:
#   if not line.startswith('#'):
#   sampleID, tumorBam, normalBam = line.strip().split()
# tumorID = sampleNameBam(tumorBam)
# normalID = sampleNameBam(normalBam)
# sample_list.append('indel_' + sampleID)
# cmd = f"indel_{sampleID} = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/strelka/{sampleID}_{tumorID}_{normalID}results/variants/somatic.indels.passonly.vep.vcf.gz.tidy.txt', sep='\\t', header=T, na.strings = ''))"
# print(cmd)
# ','.join(sample_list)
indel_4638 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/strelka/4638_4638_4361/results/variants/somatic.indels.passonly.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings = ''))
indel_5685 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/strelka/5685_5685_5497/results/variants/somatic.indels.passonly.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings = ''))
indel_3897 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/strelka/3897_3897_3665/results/variants/somatic.indels.passonly.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings = ''))
indel_4663 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/strelka/4663_4663_4385/results/variants/somatic.indels.passonly.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings = ''))
indel_4031 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/strelka/4031_4031_3801/results/variants/somatic.indels.passonly.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings = ''))
indel_4132 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/strelka/4132_4132_3904/results/variants/somatic.indels.passonly.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings = ''))
indel_4043 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/strelka/4043_4043_3813/results/variants/somatic.indels.passonly.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings = ''))
indel_3869 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/strelka/3869_3869_3637/results/variants/somatic.indels.passonly.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings = ''))
indel_4491 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/strelka/4491_4491_4261/results/variants/somatic.indels.passonly.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings = ''))
indel_3346 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/strelka/3346_3346_3122/results/variants/somatic.indels.passonly.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings = ''))
indel_4509 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/strelka/4509_4509_4279/results/variants/somatic.indels.passonly.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings = ''))
indel_5702 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/strelka/5702_5702_5487/results/variants/somatic.indels.passonly.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings = ''))
indel_4644 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/strelka/4644_4644_4367/results/variants/somatic.indels.passonly.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings = ''))
indel_4527 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/strelka/4527_4527_4298/results/variants/somatic.indels.passonly.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings = ''))
indel_4489 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/strelka/4489_4489_4259/results/variants/somatic.indels.passonly.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings = ''))
indel_4515 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/strelka/4515_4515_4285/results/variants/somatic.indels.passonly.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings = ''))
indel_2643 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/strelka/2643_2643_2409/results/variants/somatic.indels.passonly.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings = ''))
indel_4593 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/strelka/4593_4593_4454/results/variants/somatic.indels.passonly.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings = ''))
indel_3796 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/strelka/3796_3796_3570/results/variants/somatic.indels.passonly.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings = ''))
indel_2581 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/strelka/2581_2581_2574/results/variants/somatic.indels.passonly.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings = ''))
indel_5635 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/strelka/5635_5635_5486/results/variants/somatic.indels.passonly.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings = ''))
indel_2665 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/strelka/2665_2665_2433/results/variants/somatic.indels.passonly.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings = ''))
indel_2972 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/strelka/2972_2972_3007/results/variants/somatic.indels.passonly.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings = ''))
indel_3396 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/strelka/3396_3396_3175/results/variants/somatic.indels.passonly.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings = ''))
indel_3018 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/strelka/3018_3018_2913/results/variants/somatic.indels.passonly.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings = ''))
indel_2616 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/strelka/2616_2616_5349/results/variants/somatic.indels.passonly.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings = ''))
indel_3971 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/strelka/3971_3971_3740/results/variants/somatic.indels.passonly.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings = ''))
indel_5032 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/strelka/5032_5032_4747/results/variants/somatic.indels.passonly.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings = ''))
indel_3083 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/strelka/3083_3083_2857/results/variants/somatic.indels.passonly.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings = ''))
indel_2685 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/strelka/2685_2685_2459/results/variants/somatic.indels.passonly.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings = ''))
indel_3373 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/strelka/3373_3373_3151/results/variants/somatic.indels.passonly.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings = ''))
indel_2655 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/strelka/2655_2655_2422/results/variants/somatic.indels.passonly.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings = ''))
indel_4730 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/strelka/4730_4730_4451/results/variants/somatic.indels.passonly.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings = ''))
indel_5580 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/strelka/5580_5580_5488/results/variants/somatic.indels.passonly.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings = ''))
indel_4529 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/strelka/4529_4529_4300/results/variants/somatic.indels.passonly.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings = ''))
indel_5650 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/strelka/5650_5650_5439/results/variants/somatic.indels.passonly.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings = ''))
indel_2709 = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/strelka/2709_2709_2486/results/variants/somatic.indels.passonly.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings = ''))

indel_total = rbind(indel_4638,indel_5685,indel_3897,indel_4663,indel_4031,indel_4132,indel_4043,indel_3869,indel_4491,indel_3346,indel_4509,indel_5702,indel_4644,indel_4527,indel_4489,indel_4515,indel_2643,indel_4593,indel_3796,indel_2581,indel_5635,indel_2665,indel_2972,indel_3396,indel_3018,indel_2616,indel_3971,indel_5032,indel_3083,indel_2685,indel_3373,indel_2655,indel_4730,indel_5580,indel_4529,indel_5650,indel_2709)

indel_total %>% mutate(sampleName = as.character(sampleName))
indel_total = indel_total %>% mutate(sampleName = as.character(sampleName))


indel_total_count = indel_total %>% group_by(sampleName) %>% summarise(n=n())
indel_total_count$type = 'indel'
snv_indel_count = rbind(indel_total_count, combined_count)
########################################
# Plot both SNV and Indel
########################################


snv_indel_count_plot = ggplot(snv_indel_count, aes(x=sampleName, y=n, fill=factor(type, levels=c('Non-coding', 'Coding', 'indel')))) + ylab('Number of mutations') +
  geom_bar(stat = "identity") + theme_pubr()+ theme(axis.text.x=element_text(angle=90,hjust=1)) + scale_fill_discrete(name = "Mutation Types")



df_all %>% filter(str_detect(Consequence, 'missense_variant|stop_gain|stop_lost|splice'))
common_mutations = df_all %>% filter(Consequence %in% c('missense_variant', 'stop_gain', 'stop_lost', 'splice')) %>%
  group_by(sampleName, SYMBOL) %>% summarise(count =n()) %>% group_by(SYMBOL) %>% summarise(gene_occurence = n()) %>% arrange(desc(gene_occurence)) %>%
  filter(!str_detect(SYMBOL, '^IG'), gene_occurence > 2)

common_mutations_including_indel = rbind(indel_total, df_all) %>% filter(Consequence %in% c('missense_variant', 'stop_gain', 'stop_lost', 'splice')) %>%
  group_by(sampleName, SYMBOL) %>% summarise(count =n()) %>% group_by(SYMBOL) %>% summarise(gene_occurence = n()) %>% arrange(desc(gene_occurence)) %>%
  filter(!str_detect(SYMBOL, '^IG'), gene_occurence > 2)
########################################
# Landscape of commonly mutated genes
########################################
allSamples = factor(unique(df_all$sampleName))

indel_snv_all = rbind(df_all, indel_total)
landscape_data = indel_snv_all %>% filter(str_detect(Consequence, 'missense_variant|stop_gain|stop_lost|frameshift_variant')) %>%
  filter(SYMBOL %in% c( as.character(common_mutations_including_indel$SYMBOL), 'KRAS', 'FAM46C', 'BRAF', 'TP53', 'DIS3', 'PRDM1', 'SP140', 'TRAF3', 'ATM', 'CCND1', 'RB1', 'HISTH1E', 'LTB', 'IRF4', 'FGFR3', 'RB1', 'ACTG1', 'CYLD', 'MAX', 'ATR')) %>% 
  select(sampleName, SYMBOL, Consequence)
levels(landscape_data$sampleName) = allSamples


########################################
# Import SV data
########################################
sv_all = as.tibble(read.table('~/Documents/cjyoon/Projects/myeloma/analysis/sv_intersect/mm37_all_filtered_sv.tsv', col.names = c('bp1', 'bp2', 'sv_type', 'sample')))
sv_count = sv_all %>% separate(sample, c('sampleID', 'dummy'), extra='drop') %>% group_by(sampleID, sv_type) %>% tally() 
add_row(sv_count, sampleID='2616', sv_type='BND', n=0)
# 2616 has no overlapping SV in manta and delly. thus empty file. but need to show this in plot
sv_count_0sv = tribble(~sampleID, ~sv_type, ~n, 
                       "2616", "BND", 0, 
                       "2616", "DUP", 0, 
                       "2616", "INV", 0, 
                       "2616", "DEL", 0)
sv_count_0sv$n = as.integer(sv_count_0sv$n)
sv_count_all = bind_rows(sv_count, sv_count_0sv)

sv_count_plot = sv_count_all%>% 
  ggplot(aes(x=sampleID, y=n, fill=sv_type))+ geom_bar(stat='identity') + theme_pubr() + 
  theme(axis.text.x=element_text(angle=90, hjust=1)) +scale_fill_discrete(name = "SV Types") +
  ylab('Number of mutations') 
# median SV event
sv_count_bysample = (sv_all %>% separate(sample, c('sampleID', 'dummy'), extra='drop') %>% group_by(sampleID, sv_type) %>% group_by(sampleID) %>% tally())
median(sv_count_bysample$n) # median SV event per sample = 30 



snv_landscape = landscape_data %>% 
  ggplot(aes(x=sampleName, y=SYMBOL, fill=Consequence)) + geom_tile() + theme_pubr() + theme(axis.text.x=element_text(angle=90, hjust=1)) + xlab('SampleID') + ylab('Gene') 

# these sampleIDs do not have SNVs in the selected genes. BUT to align with panel A and B with SNV/SV counts will manually add as NA values
sampleName = c("2581", "2616", "2972", "3018", "3396", "4031", "4132", "4644", "4663", "5032", 'dummy') # dummy was added to prevent column type from becoming 'logical' value which prevent 'rbind'
SYMBOL = c(rep(NA, 10), 'bla')
Consequence = c(rep(NA, 10), 'bla')


missing_from_landscape = as.tibble(data.frame(sampleName, SYMBOL, Consequence, stringsAsFactors=FALSE))
# Combine landcape data with samples without mutation in these genes 

landscape_all = rbind(missing_from_landscape, landscape_data)
snv_landscape_all = landscape_all %>% filter(sampleName != 'dummy') %>% 
  ggplot(aes(x=sampleName, y=SYMBOL, fill=Consequence)) + geom_tile(color='black', size=0.5) + theme_pubr() + theme(axis.text.x=element_text(angle=90, hjust=1)) + xlab('SampleID') + ylab('Gene')


########################################
# Combine counts and Landscape into a single plot
########################################
sv_snv_plot = plot_grid(snv_indel_count_plot, sv_count_plot, snv_landscape_all, labels = c('A', 'B', 'C'), ncol=1, align='v', rel_heights = c(1, 1, 2))
ggsave('~/Documents/cjyoon/Projects/myeloma/figures/landscape_w_indel_v3_2018.07.15.png', width=10, height=20, dpi=200)


