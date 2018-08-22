# Script to see if the distribution of HLA class I alleles are different from Korean population
# Korean reference HLA allele frequency was obtained from 613 healthy Koreans.
# In et al (2015) Ann Lab Medicine 10.3343/alm.2015.35.4.429
# 2018.07.12 CJY
setwd("~/Documents/cjyoon/Projects/myeloma/analysis/polysolver/analysis")

library(tidyverse)
library(ggpubr)

n = 37 * 2 # total number of people * 2 alleles per person


##############################################################################################################
# HLA-A
hla_a_observed = read_delim('hla_a_observed_mm.txt', delim='\t')
hla_a_korean = as.tibble(read.table('hla_a_korean.txt', header=T, stringsAsFactors=FALSE))
hla_a_expected = hla_a_korean %>% mutate(expected = allele_frequency / 100 * n)
hla_a_table = hla_a_expected %>% full_join(hla_a_observed) %>% filter(observed != 'NA' & expected != 'NA')
chisq.test(x=hla_a_table$observed, y=hla_a_table$expected)

# Pearson's Chi-squared test
# 
# data:  hla_a_table$observed and hla_a_table$expected
# X-squared = 105, df = 98, p-value = 0.296

hla_a_table_observed = hla_a_table %>% gather(key=type, value=count, c('expected', 'observed')) %>% filter(type=='observed')
hla_a_table_expected = hla_a_table %>% gather(key=type, value=count, c('expected', 'observed')) %>% filter(type=='expected')

# Draw Expected and Observed HLA allele frequencies
hla_a_table_observed %>% 
  ggplot(aes(x=reorder(HLA_allele, -count), y=count)) + 
  geom_bar(stat="identity") + 
  geom_point(data=hla_a_table_expected, aes(x=HLA_allele, y=count)) +
  theme_pubr() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  labs(x='HLA-A alleles', y='Allele counts')+ 
  annotate(geom="text", x=Inf, y=Inf, hjust =1, vjust=1, label="Chi-square\np-value=0.296",color="red")


  
##############################################################################################################
# HLA-B
hla_b_observed = read_delim('hla_b_observed_mm.txt', delim='\t')
hla_b_korean = as.tibble(read.table('hla_b_korean.txt', header=T, stringsAsFactors=FALSE))
hla_b_expected = hla_b_korean %>% mutate(expected = allele_frequency / 100 * n)
hla_b_table = hla_b_expected %>% full_join(hla_b_observed)  %>% filter(observed != 'NA' & expected != 'NA')
chisq.test(x=hla_b_table$observed, y=hla_b_table$expected)
# Pearson's Chi-squared test
# 
# data:  hla_b_table$observed and hla_b_table$expected
# X-squared = 171.11, df = 147, p-value = 0.08467

hla_b_table_observed = hla_b_table %>% gather(key=type, value=count, c('expected', 'observed')) %>% filter(type=='observed')
hla_b_table_expected = hla_b_table %>% gather(key=type, value=count, c('expected', 'observed')) %>% filter(type=='expected')

# Draw Expected and Observed HLA allele frequencies
hla_b_table_observed %>% 
  ggplot(aes(x=reorder(HLA_allele, -count), y=count)) + 
  geom_bar(stat="identity") + 
  geom_point(data=hla_b_table_expected, aes(x=HLA_allele, y=count)) +
  theme_pubr() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  labs(x='HLA-B alleles', y='Allele counts')+ 
  annotate(geom="text", x=Inf, y=Inf, hjust =1, vjust=1, label="Chi-square\np-value=0.08467",color="red")


##############################################################################################################
# HLA-C
hla_c_observed = read_delim('hla_c_observed_mm.txt', delim='\t')
hla_c_korean = as.tibble(read.table('hla_c_korean.txt', header=T, stringsAsFactors=FALSE))
hla_c_expected = hla_c_korean %>% mutate(expected = allele_frequency / 100 * n)
hla_c_table = hla_c_expected %>% full_join(hla_c_observed) %>% filter(observed != 'NA' & expected != 'NA')
chisq.test(x=hla_c_table$observed, y=hla_c_table$expected)
# Pearson's Chi-squared test
# 
# data:  hla_c_table$observed and hla_c_table$expected
# X-squared = 144, df = 136, p-value = 0.3029

hla_c_table_observed = hla_c_table %>% gather(key=type, value=count, c('expected', 'observed')) %>% filter(type=='observed')
hla_c_table_expected = hla_c_table %>% gather(key=type, value=count, c('expected', 'observed')) %>% filter(type=='expected')

# Draw Expected and Observed HLA allele frequencies
hla_c_table_observed %>% 
  ggplot(aes(x=reorder(HLA_allele, -count), y=count)) + 
  geom_bar(stat="identity") + 
  geom_point(data=hla_c_table_expected, aes(x=HLA_allele, y=count)) +
  theme_pubr() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  labs(x='HLA-C alleles', y='Allele counts') + 
  annotate(geom="text", x=Inf, y=Inf, hjust =1, vjust=1, label="Chi-square\np-value=0.3029",color="red")

