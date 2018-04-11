# Summarize SNVs in Multiple Myeloma WGS (n=29)
# Jongsoo Yoon (cjyoon@kaist.ac.kr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggbio)

########################################
# Import and Combine tidied data
########################################

#python code to generate those lines 
# for i in [re.sub('_mns.vep.vcf.gz.tidy.txt', '', sampleID) for sampleID in os.listdir('/home/users/cjyoon/Projects/myeloma/analysis/annotated_snvs/')]:
# fn = f"df_{i} = as.tibble(read.table('~/Projects/myeloma/analysis/annotated_snvs/{i}_mns.vep.vcf.gz.tidy.txt', sep='\t', header=T, na.strings=''))"
# print(fn)
# fn = f"df_{i}$sampleID = '{i}'"
# print(fn)

df_2665 = as.tibble(read.table('~/Projects/myeloma/analysis/annotated_snvs/2665_mns.vep.vcf.gz.tidy.txt', sep='	', header=T, na.strings=""))
df_2665$sampleID = '2665'
df_3018 = as.tibble(read.table('~/Projects/myeloma/analysis/annotated_snvs/3018_mns.vep.vcf.gz.tidy.txt', sep='	', header=T, na.strings=""))
df_3018$sampleID = '3018'
df_4489 = as.tibble(read.table('~/Projects/myeloma/analysis/annotated_snvs/4489_mns.vep.vcf.gz.tidy.txt', sep='	', header=T, na.strings=""))
df_4489$sampleID = '4489'
df_4509 = as.tibble(read.table('~/Projects/myeloma/analysis/annotated_snvs/4509_mns.vep.vcf.gz.tidy.txt', sep='	', header=T, na.strings=""))
df_4509$sampleID = '4509'
df_3396 = as.tibble(read.table('~/Projects/myeloma/analysis/annotated_snvs/3396_mns.vep.vcf.gz.tidy.txt', sep='	', header=T, na.strings=""))
df_3396$sampleID = '3396'
df_4491 = as.tibble(read.table('~/Projects/myeloma/analysis/annotated_snvs/4491_mns.vep.vcf.gz.tidy.txt', sep='	', header=T, na.strings=""))
df_4491$sampleID = '4491'
df_4031 = as.tibble(read.table('~/Projects/myeloma/analysis/annotated_snvs/4031_mns.vep.vcf.gz.tidy.txt', sep='	', header=T, na.strings=""))
df_4031$sampleID = '4031'
df_4644 = as.tibble(read.table('~/Projects/myeloma/analysis/annotated_snvs/4644_mns.vep.vcf.gz.tidy.txt', sep='	', header=T, na.strings=""))
df_4644$sampleID = '4644'
df_4132 = as.tibble(read.table('~/Projects/myeloma/analysis/annotated_snvs/4132_mns.vep.vcf.gz.tidy.txt', sep='	', header=T, na.strings=""))
df_4132$sampleID = '4132'
df_4527 = as.tibble(read.table('~/Projects/myeloma/analysis/annotated_snvs/4527_mns.vep.vcf.gz.tidy.txt', sep='	', header=T, na.strings=""))
df_4527$sampleID = '4527'
df_3971 = as.tibble(read.table('~/Projects/myeloma/analysis/annotated_snvs/3971_mns.vep.vcf.gz.tidy.txt', sep='	', header=T, na.strings=""))
df_3971$sampleID = '3971'
df_2709 = as.tibble(read.table('~/Projects/myeloma/analysis/annotated_snvs/2709_mns.vep.vcf.gz.tidy.txt', sep='	', header=T, na.strings=""))
df_2709$sampleID = '2709'
df_4043 = as.tibble(read.table('~/Projects/myeloma/analysis/annotated_snvs/4043_mns.vep.vcf.gz.tidy.txt', sep='	', header=T, na.strings=""))
df_4043$sampleID = '4043'
df_3373 = as.tibble(read.table('~/Projects/myeloma/analysis/annotated_snvs/3373_mns.vep.vcf.gz.tidy.txt', sep='	', header=T, na.strings=""))
df_3373$sampleID = '3373'
df_2655 = as.tibble(read.table('~/Projects/myeloma/analysis/annotated_snvs/2655_mns.vep.vcf.gz.tidy.txt', sep='	', header=T, na.strings=""))
df_2655$sampleID = '2655'
df_5032 = as.tibble(read.table('~/Projects/myeloma/analysis/annotated_snvs/5032_mns.vep.vcf.gz.tidy.txt', sep='	', header=T, na.strings=""))
df_5032$sampleID = '5032'
df_2643 = as.tibble(read.table('~/Projects/myeloma/analysis/annotated_snvs/2643_mns.vep.vcf.gz.tidy.txt', sep='	', header=T, na.strings=""))
df_2643$sampleID = '2643'
df_3346 = as.tibble(read.table('~/Projects/myeloma/analysis/annotated_snvs/3346_mns.vep.vcf.gz.tidy.txt', sep='	', header=T, na.strings=""))
df_3346$sampleID = '3346'
df_4730 = as.tibble(read.table('~/Projects/myeloma/analysis/annotated_snvs/4730_mns.vep.vcf.gz.tidy.txt', sep='	', header=T, na.strings=""))
df_4730$sampleID = '4730'
df_3897 = as.tibble(read.table('~/Projects/myeloma/analysis/annotated_snvs/3897_mns.vep.vcf.gz.tidy.txt', sep='	', header=T, na.strings=""))
df_3897$sampleID = '3897'
df_3796 = as.tibble(read.table('~/Projects/myeloma/analysis/annotated_snvs/3796_mns.vep.vcf.gz.tidy.txt', sep='	', header=T, na.strings=""))
df_3796$sampleID = '3796'
df_4593 = as.tibble(read.table('~/Projects/myeloma/analysis/annotated_snvs/4593_mns.vep.vcf.gz.tidy.txt', sep='	', header=T, na.strings=""))
df_4593$sampleID = '4593'
df_2685 = as.tibble(read.table('~/Projects/myeloma/analysis/annotated_snvs/2685_mns.vep.vcf.gz.tidy.txt', sep='	', header=T, na.strings=""))
df_2685$sampleID = '2685'
df_3869 = as.tibble(read.table('~/Projects/myeloma/analysis/annotated_snvs/3869_mns.vep.vcf.gz.tidy.txt', sep='	', header=T, na.strings=""))
df_3869$sampleID = '3869'
df_4663 = as.tibble(read.table('~/Projects/myeloma/analysis/annotated_snvs/4663_mns.vep.vcf.gz.tidy.txt', sep='	', header=T, na.strings=""))
df_4663$sampleID = '4663'
df_4529 = as.tibble(read.table('~/Projects/myeloma/analysis/annotated_snvs/4529_mns.vep.vcf.gz.tidy.txt', sep='	', header=T, na.strings=""))
df_4529$sampleID = '4529'
df_4638 = as.tibble(read.table('~/Projects/myeloma/analysis/annotated_snvs/4638_mns.vep.vcf.gz.tidy.txt', sep='	', header=T, na.strings=""))
df_4638$sampleID = '4638'
df_4515 = as.tibble(read.table('~/Projects/myeloma/analysis/annotated_snvs/4515_mns.vep.vcf.gz.tidy.txt', sep='	', header=T, na.strings=""))
df_4515$sampleID = '4515'
df_3083 = as.tibble(read.table('~/Projects/myeloma/analysis/annotated_snvs/3083_mns.vep.vcf.gz.tidy.txt', sep='	', header=T, na.strings=""))
df_3083$sampleID = '3083'

# Combined DataFrame for all SNVs
df_all = rbind(df_2665,df_3018,df_4489,df_4509,df_3396,df_4491,df_4031,df_4644,df_4132,df_4527,df_3971,df_2709,df_4043,df_3373,df_2655,df_5032,df_2643,df_3346,df_4730,df_3897,df_3796,df_4593,df_2685,df_3869,df_4663,df_4529,df_4638,df_4515,df_3083)


########################################
# Number of SNVs in Exons and WGS
########################################

noncoding_mutation_count = df_all %>% select(SYMBOL, sampleID) %>% filter(is.na(SYMBOL)) %>% group_by(sampleID) %>% tally()
mutation_count$type='Non-coding'
coding_mutation_count = df_all %>% select(SYMBOL, sampleID) %>% filter(!(is.na(SYMBOL))) %>% group_by(sampleID) %>% tally() 
coding_mutation_count$type='Coding'

combined_count = rbind(mutation_count, coding_mutation_count)
combined_count$sampleID = factor(combined_count$sampleID, levels = rev(levels(combined_count$sampleID)))
ggplot(combined_count, aes(x=sampleID, y=n, fill=factor(type, levels=c('Non-coding', 'Coding')))) + ylab('Number of mutations') +
         geom_bar(stat = "identity") + theme_pubr()+ theme(axis.text.x=element_text(angle=90,hjust=1)) + scale_fill_discrete(name = "Mutation Types")



########################################
# Import SV data
########################################
sv_all = as.tibble(read.table('/home/users/cjyoon/Projects/myeloma/analysis/sv_intersect/sv_filtered_29mm_all.tsv', col.names = c('bp2', 'bp2', 'sv_type', 'sample')))
sv_all %>% separate(sample, c('sampleID', 'dummy'), extra='drop')
