# compares SNU's cibersort result vs Tumor cell fraction results from sequenza
# 2018.07.29 CJY

library(tidyverse)
library(readxl)
library(ggplot2)
cibersort_path = '~/Documents/cjyoon/Projects/myeloma/analysis/snu_data/CIBERSORT.output_paired.MM_n82.txt'
cibersort_df = read_delim(cibersort_path, delim='\t')
cibersort_df
cibersort_plasma = cibersort_df %>% separate(`Input Sample`, into=c('sample', 'sort_status')) %>% filter(sort_status == 'positive') %>% dplyr::select(c('sample', 'Plasma cells', 'B cells naive', 'B cells memory'))


sample_summary_path = '~/Documents/julab/projects/myeloma/sample_summary.xlsx'
sample_summary_df = read_excel(sample_summary_path)
sample_summary_df$Tumor = as.character(sample_summary_df$Tumor)

compare_cibersort_sequenza = sample_summary_df %>% dplyr::select(c('Tumor', 'cellularity')) %>% left_join(cibersort_plasma, by=c('Tumor'='sample'))
colnames(compare_cibersort_sequenza) = c('sampleID', 'sequenza_cellularity', 'cibersort_plasma_cells', 'cibersort_naiveB', 'cibersort_memoryB')
compare_cibersort_sequenza = compare_cibersort_sequenza %>% mutate(cibersort_allB = cibersort_plasma_cells + cibersort_naiveB + cibersort_memoryB)

ggplot(compare_cibersort_sequenza, aes(x=sequenza_cellularity, y=cibersort_plasma_cells)) + geom_point() +
  geom_smooth(method='lm') + xlim(c(0, 1)) + ylim(c(0, 1))
ggsave('~/Documents/cjyoon/Projects/myeloma/analysis/snu_data/sequenza_cellularity_vs_cibersort_plasmacell.png')

ggplot(compare_cibersort_sequenza, aes(x=sequenza_cellularity, y=cibersort_naiveB)) + geom_point() +
  geom_smooth(method='lm') + xlim(c(0, 1)) + ylim(c(0, 1))
ggsave('~/Documents/cjyoon/Projects/myeloma/analysis/snu_data/sequenza_cellularity_vs_cibersort_naiveB.png')

ggplot(compare_cibersort_sequenza, aes(x=sequenza_cellularity, y=cibersort_memoryB)) + geom_point() +
  geom_smooth(method='lm') + xlim(c(0, 1)) + ylim(c(0, 1))
ggsave('~/Documents/cjyoon/Projects/myeloma/analysis/snu_data/sequenza_cellularity_vs_cibersort_memoryB.png')


ggplot(compare_cibersort_sequenza, aes(x=sequenza_cellularity, y=cibersort_allB)) + geom_point() +
  geom_smooth(method='lm') + xlim(c(0, 1)) + ylim(c(0, 1))
ggsave('~/Documents/cjyoon/Projects/myeloma/analysis/snu_data/sequenza_cellularity_vs_cibersort_allB.png')

