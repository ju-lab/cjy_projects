# Summarize SNVs in Multiple Myeloma WGS (n=29)
# Jongsoo Yoon (cjyoon@kaist.ac.kr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggbio)
library(cowplot)
# Add Anaconda installed Library paths
.libPaths( c( .libPaths(), "/home/users/cjyoon/anaconda3/lib/R/library") )
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

# Median value 
snv_count_bysample = df_all %>% group_by(sampleID) %>% tally() 
median(snv_count_bysample$n) # 6302

snv_count = ggplot(combined_count, aes(x=sampleID, y=n, fill=factor(type, levels=c('Non-coding', 'Coding')))) + ylab('Number of mutations') +
         geom_bar(stat = "identity") + theme_pubr()+ theme(axis.text.x=element_text(angle=90,hjust=1)) + scale_fill_discrete(name = "Mutation Types")



########################################
# Import SV data
########################################
sv_all = as.tibble(read.table('/home/users/cjyoon/Projects/myeloma/analysis/sv_intersect/sv_filtered_29mm_all.tsv', col.names = c('bp1', 'bp2', 'sv_type', 'sample')))
sv_count = sv_all %>% separate(sample, c('sampleID', 'dummy'), extra='drop') %>% group_by(sampleID, sv_type) %>% tally() %>% 
  ggplot(aes(x=sampleID, y=n, fill=sv_type))+ geom_bar(stat='identity') + theme_pubr() + 
  theme(axis.text.x=element_text(angle=90, hjust=1)) +scale_fill_discrete(name = "SV Types") +
  ylab('Number of mutations')

# median SV event
sv_count_bysample = (sv_all %>% separate(sample, c('sampleID', 'dummy'), extra='drop') %>% group_by(sampleID, sv_type) %>% group_by(sampleID) %>% tally())
median(sv_count_bysample$n) # median SV event per sample = 28

########################################
# Combine counts into a single plot
########################################
sv_snv_plot = plot_grid(snv_count, sv_count, labels = c('A', 'B'), ncol=1, align='v')
ggsave(filename = '/home/users/cjyoon/Projects/myeloma/figures/sv_snv_counts.png', plot = sv_snv_plot, dpi = 200)



########################################
# Import Clinical Data
########################################
clinical_df = as.tibble(read.table('/home/users/cjyoon/Projects/myeloma/sample_info/clinical_data.txt', header=T, colClasses=c('character', 'character', 'integer', 'factor', 'double', 'character')))

########################################
# Dependence on Age
# -> Negative 
########################################
my.formula <- y ~ x
sv_count_bysample %>% left_join(clinical_df, by = c('sampleID'='sampleName')) %>% 
  ggplot(aes(x=Age, y=n)) + geom_point() + geom_smooth(method='lm',formula=y~x) +  
  stat_poly_eq(formula = my.formula, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),  parse = TRUE)
snv_count_bysample %>% left_join(clinical_df, by = c('sampleID'='sampleName')) %>% 
  ggplot(aes(x=Age, y=n)) + geom_point() + geom_smooth(method='lm',formula=y~x) +  
  stat_poly_eq(formula = my.formula, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),  parse = TRUE)


# Correlation between SNV and SV count 
sv_count_by_sample_join = sv_count_bysample %>% left_join(clinical_df, by = c('sampleID'='sampleName')) %>% mutate(sv_count = n) %>% select(sampleID, sv_count)
snv_count_by_sample_join = snv_count_bysample %>% left_join(clinical_df, by = c('sampleID'='sampleName')) %>% mutate(snv_count = n) %>% select(sampleID, snv_count)

sv_snv_correlation = sv_count_by_sample_join %>% left_join(snv_count_by_sample_join) %>%
  ggplot(aes(x=snv_count, y=sv_count)) + geom_point() + geom_smooth(method='lm', formulat=y~x) + 
  stat_poly_eq(formula = my.formula, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),  parse = TRUE) + 
  xlab('SNV count') + ylab('SV count') 
ggsave(filename='/home/users/cjyoon/Projects/myeloma/figures/sv_snv_correlation.png',sv_snv_correlation , dpi=200 )

########################################
# Landscape plot 
########################################


df_all %>% filter(Consequence == 'intron_variant') %>% group_by(SYMBOL) %>% tally() %>% arrange(desc(n))
# # A tibble: 11,170 x 2
# SYMBOL      n
# <chr>   <int>
#   1 ROBO2     417
# 2 CSMD1     400
# 3 CNTNAP2   371
# 4 EYS       367
# 5 LRP1B     364
# 6 PTPRD     324
# 7 CTNNA2    303
# 8 CSMD3     289
# 9 GRID2     286
# 10 CDH12     272

df_all %>% filter(Consequence != 'intron_variant') %>% group_by(SYMBOL) %>% tally() %>% arrange(desc(n))
