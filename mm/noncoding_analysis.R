library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggbio)
library(cowplot)
library(stringr)

df_2665 = as.tibble(read.table('~/Projects/myeloma/analysis/annotated_snvs/withvaf/2665_mns.vep.vcf.gz.tidy.txt', sep='	', header=T, na.strings=''))
df_3018 = as.tibble(read.table('~/Projects/myeloma/analysis/annotated_snvs/withvaf/3018_mns.vep.vcf.gz.tidy.txt', sep='	', header=T, na.strings=''))
df_4489 = as.tibble(read.table('~/Projects/myeloma/analysis/annotated_snvs/withvaf/4489_mns.vep.vcf.gz.tidy.txt', sep='	', header=T, na.strings=''))
df_4509 = as.tibble(read.table('~/Projects/myeloma/analysis/annotated_snvs/withvaf/4509_mns.vep.vcf.gz.tidy.txt', sep='	', header=T, na.strings=''))
df_3396 = as.tibble(read.table('~/Projects/myeloma/analysis/annotated_snvs/withvaf/3396_mns.vep.vcf.gz.tidy.txt', sep='	', header=T, na.strings=''))
df_4491 = as.tibble(read.table('~/Projects/myeloma/analysis/annotated_snvs/withvaf/4491_mns.vep.vcf.gz.tidy.txt', sep='	', header=T, na.strings=''))
df_4031 = as.tibble(read.table('~/Projects/myeloma/analysis/annotated_snvs/withvaf/4031_mns.vep.vcf.gz.tidy.txt', sep='	', header=T, na.strings=''))
df_4644 = as.tibble(read.table('~/Projects/myeloma/analysis/annotated_snvs/withvaf/4644_mns.vep.vcf.gz.tidy.txt', sep='	', header=T, na.strings=''))
df_4132 = as.tibble(read.table('~/Projects/myeloma/analysis/annotated_snvs/withvaf/4132_mns.vep.vcf.gz.tidy.txt', sep='	', header=T, na.strings=''))
df_4527 = as.tibble(read.table('~/Projects/myeloma/analysis/annotated_snvs/withvaf/4527_mns.vep.vcf.gz.tidy.txt', sep='	', header=T, na.strings=''))
df_3971 = as.tibble(read.table('~/Projects/myeloma/analysis/annotated_snvs/withvaf/3971_mns.vep.vcf.gz.tidy.txt', sep='	', header=T, na.strings=''))
df_2709 = as.tibble(read.table('~/Projects/myeloma/analysis/annotated_snvs/withvaf/2709_mns.vep.vcf.gz.tidy.txt', sep='	', header=T, na.strings=''))
df_4043 = as.tibble(read.table('~/Projects/myeloma/analysis/annotated_snvs/withvaf/4043_mns.vep.vcf.gz.tidy.txt', sep='	', header=T, na.strings=''))
df_3373 = as.tibble(read.table('~/Projects/myeloma/analysis/annotated_snvs/withvaf/3373_mns.vep.vcf.gz.tidy.txt', sep='	', header=T, na.strings=''))
df_2655 = as.tibble(read.table('~/Projects/myeloma/analysis/annotated_snvs/withvaf/2655_mns.vep.vcf.gz.tidy.txt', sep='	', header=T, na.strings=''))
df_5032 = as.tibble(read.table('~/Projects/myeloma/analysis/annotated_snvs/withvaf/5032_mns.vep.vcf.gz.tidy.txt', sep='	', header=T, na.strings=''))
df_2643 = as.tibble(read.table('~/Projects/myeloma/analysis/annotated_snvs/withvaf/2643_mns.vep.vcf.gz.tidy.txt', sep='	', header=T, na.strings=''))
df_3346 = as.tibble(read.table('~/Projects/myeloma/analysis/annotated_snvs/withvaf/3346_mns.vep.vcf.gz.tidy.txt', sep='	', header=T, na.strings=''))
df_4730 = as.tibble(read.table('~/Projects/myeloma/analysis/annotated_snvs/withvaf/4730_mns.vep.vcf.gz.tidy.txt', sep='	', header=T, na.strings=''))
df_3897 = as.tibble(read.table('~/Projects/myeloma/analysis/annotated_snvs/withvaf/3897_mns.vep.vcf.gz.tidy.txt', sep='	', header=T, na.strings=''))
df_3796 = as.tibble(read.table('~/Projects/myeloma/analysis/annotated_snvs/withvaf/3796_mns.vep.vcf.gz.tidy.txt', sep='	', header=T, na.strings=''))
df_4593 = as.tibble(read.table('~/Projects/myeloma/analysis/annotated_snvs/withvaf/4593_mns.vep.vcf.gz.tidy.txt', sep='	', header=T, na.strings=''))
df_2685 = as.tibble(read.table('~/Projects/myeloma/analysis/annotated_snvs/withvaf/2685_mns.vep.vcf.gz.tidy.txt', sep='	', header=T, na.strings=''))
df_3869 = as.tibble(read.table('~/Projects/myeloma/analysis/annotated_snvs/withvaf/3869_mns.vep.vcf.gz.tidy.txt', sep='	', header=T, na.strings=''))
df_4663 = as.tibble(read.table('~/Projects/myeloma/analysis/annotated_snvs/withvaf/4663_mns.vep.vcf.gz.tidy.txt', sep='	', header=T, na.strings=''))
df_4529 = as.tibble(read.table('~/Projects/myeloma/analysis/annotated_snvs/withvaf/4529_mns.vep.vcf.gz.tidy.txt', sep='	', header=T, na.strings=''))
df_4638 = as.tibble(read.table('~/Projects/myeloma/analysis/annotated_snvs/withvaf/4638_mns.vep.vcf.gz.tidy.txt', sep='	', header=T, na.strings=''))
df_4515 = as.tibble(read.table('~/Projects/myeloma/analysis/annotated_snvs/withvaf/4515_mns.vep.vcf.gz.tidy.txt', sep='	', header=T, na.strings=''))
df_3083 = as.tibble(read.table('~/Projects/myeloma/analysis/annotated_snvs/withvaf/3083_mns.vep.vcf.gz.tidy.txt', sep='	', header=T, na.strings=''))

df_all = rbind(df_2665,df_3018,df_4489,df_4509,df_3396,df_4491,df_4031,df_4644,df_4132,df_4527,df_3971,df_2709,df_4043,df_3373,df_2655,df_5032,df_2643,df_3346,df_4730,df_3897,df_3796,df_4593,df_2685,df_3869,df_4663,df_4529,df_4638,df_4515,df_3083)
df_all = df_all %>% mutate(sampleName = as.character(sampleName))

# per each chromosome, see if any region show increased density of noncoding mutations
for(chr in c(1:22, 'X', 'Y')){
  print(chr)
  chr_mutations = df_all %>% filter(chromosome==as.character(chr))
  chr_noncoding = chr_mutations %>% filter(grepl('upstream_gene_variant|downstream_gene_variant|TF_binding_site_variant|regulatory_region_variant|intergenic|3_prime_UTR_variant|5_prime_UTR_variant', Consequence))  
  ggplot(chr_noncoding, aes(chr_noncoding$position)) + 
    geom_histogram(breaks=seq(0, max(chr_noncoding$position), by = 10000)) + 
    geom_density(col=2) + 
    ggtitle(str_c('Noncoding mutation densities in chr', chr)) + 
    xlab('Position')
  ggsave(str_c('/home/users/cjyoon/Projects/myeloma/analysis/noncoding/mm_noncoding_chr', as.character(chr), '.density.png'))
}

