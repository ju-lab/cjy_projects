# This script is to compare Yonsei's RNA seq result with Brainspan RNA-seq data
# Jongsoo Yoon (cjyoon@kaist.ac.kr) April 4 2018

library(ggplot2)
library(tidyverse)
############################################################################################

# Brainspan data and metadata
# Brainspan data in RPKM value
brainspan_rowmeta = as.tibble(read.table('~/Projects/yonsei_organoid/analysis/rows_metadata.csv', sep=',', header=T))
brainspan_columnmeta = as.tibble(read.table('~/Projects/yonsei_organoid/analysis/columns_metadata.csv', sep=',', header=T))    
expression = read_csv("~/Projects/yonsei_organoid/analysis/expression_matrix.csv", col_names = c('index', 1:524))

expression_w_ensembl = expression %>% left_join(brainspan_rowmeta, by=c('index'='row_num'))

############################################################################################
# Now import RSEM data for each conditions
# RSEM data in FPKM value, but since we're comparing the ordrers only, it doesn't matter
# Data read and parsing script for RSEM
format_rsem_ensembl_gene <- function(rsem_output){
  df = read_delim(rsem_output, delim='\t')
  df %>% separate(gene_id, into=c('ensembl_gene_id', 'ensembl_gene_id_version'))
}
############################################################################################
# BEM 


bem_1 = format_rsem_ensembl_gene('~/Projects/yonsei_organoid/rsem/BEM_1_rsem.genes.results')
bem_2 = format_rsem_ensembl_gene('~/Projects/yonsei_organoid/rsem/BEM_2_rsem.genes.results')
bem_3 = format_rsem_ensembl_gene('~/Projects/yonsei_organoid/rsem/BEM_3_rsem.genes.results')


bem_1$sample_name = 'bem_1'
bem_2$sample_name = 'bem_2'
bem_3$sample_name = 'bem_3'

# distinct removes duplicate rows for each rsem results (one example of duplicate row had ensembl_gene_id=IGH)
bem_1_data = bem_1 %>% select(ensembl_gene_id, FPKM, sample_name) %>% distinct(ensembl_gene_id, FPKM, sample_name)
bem_2_data = bem_2 %>% select(ensembl_gene_id, FPKM, sample_name) %>% distinct(ensembl_gene_id, FPKM, sample_name)
bem_3_data = bem_3 %>% select(ensembl_gene_id, FPKM, sample_name) %>% distinct(ensembl_gene_id, FPKM, sample_name)
bem_merged = rbind(bem_1_data, bem_2_data, bem_3_data)

bem_fpkm = bem_merged %>% group_by(ensembl_gene_id) %>% summarise(bem = mean(FPKM))
############################################################################################
# MAT
mat_1 = format_rsem_ensembl_gene('~/Projects/yonsei_organoid/rsem/MAT_1_rsem.genes.results')
mat_2 = format_rsem_ensembl_gene('~/Projects/yonsei_organoid/rsem/MAT_2_rsem.genes.results')
mat_3 = format_rsem_ensembl_gene('~/Projects/yonsei_organoid/rsem/MAT_3_rsem.genes.results')

mat_1$sample_name = 'mat_1'
mat_2$sample_name = 'mat_2'
mat_3$sample_name = 'mat_3'
mat_1_data = mat_1 %>% select(ensembl_gene_id, FPKM, sample_name) %>% distinct(ensembl_gene_id, FPKM, sample_name)
mat_2_data = mat_2 %>% select(ensembl_gene_id, FPKM, sample_name) %>% distinct(ensembl_gene_id, FPKM, sample_name)
mat_3_data = mat_3 %>% select(ensembl_gene_id, FPKM, sample_name) %>% distinct(ensembl_gene_id, FPKM, sample_name)

mat_merged = rbind(mat_1_data, mat_2_data, mat_3_data)
mat_fpkm = mat_merged %>% group_by(ensembl_gene_id) %>% summarise(mat = mean(FPKM))
############################################################################################
# NSC
nsc_1 = format_rsem_ensembl_gene('~/Projects/yonsei_organoid/rsem/NSC_1_rsem.genes.results')
nsc_2 = format_rsem_ensembl_gene('~/Projects/yonsei_organoid/rsem/NSC_2_rsem.genes.results')
nsc_3 = format_rsem_ensembl_gene('~/Projects/yonsei_organoid/rsem/NSC_3_rsem.genes.results')

nsc_1$sample_name = 'nsc_1'
nsc_2$sample_name = 'nsc_2'
nsc_3$sample_name = 'nsc_3'

nsc_1_data = nsc_1 %>% select(ensembl_gene_id, FPKM, sample_name) %>% distinct(ensembl_gene_id, FPKM, sample_name)
nsc_2_data = nsc_2 %>% select(ensembl_gene_id, FPKM, sample_name) %>% distinct(ensembl_gene_id, FPKM, sample_name)
nsc_3_data = nsc_3 %>% select(ensembl_gene_id, FPKM, sample_name) %>% distinct(ensembl_gene_id, FPKM, sample_name)

nsc_merged = rbind(nsc_1_data, nsc_2_data, nsc_3_data)
nsc_fpkm = nsc_merged %>% group_by(ensembl_gene_id) %>% summarise(nsc = mean(FPKM))
############################################################################################
# Tissue

tissue_1 = format_rsem_ensembl_gene('~/Projects/yonsei_organoid/rsem/Tissue_1_rsem.genes.results')
tissue_2 = format_rsem_ensembl_gene('~/Projects/yonsei_organoid/rsem/Tissue_2_rsem.genes.results')
tissue_3 = format_rsem_ensembl_gene('~/Projects/yonsei_organoid/rsem/Tissue_3_rsem.genes.results')

tissue_1$sample_name = 'tissue_1'
tissue_2$sample_name = 'tissue_2'
tissue_3$sample_name = 'tissue_3'

tissue_1_data = tissue_1 %>% select(ensembl_gene_id, FPKM, sample_name) %>% distinct(ensembl_gene_id, FPKM, sample_name)
tissue_2_data = tissue_2 %>% select(ensembl_gene_id, FPKM, sample_name) %>% distinct(ensembl_gene_id, FPKM, sample_name)
tissue_3_data = tissue_3 %>% select(ensembl_gene_id, FPKM, sample_name) %>% distinct(ensembl_gene_id, FPKM, sample_name)

tissue_merged = rbind(tissue_1_data, tissue_2_data, tissue_3_data)
tissue_fpkm = tissue_merged %>% group_by(ensembl_gene_id) %>% summarise(tissue = mean(FPKM))
############################################################################################
# Now combine into a single tibble for correlation calculation

expression_with_bem_mat = expression_w_ensembl %>% left_join(mat_fpkm) %>% left_join(bem_fpkm) %>% left_join(nsc_fpkm) %>% left_join(tissue_fpkm)

#Calculate pairwise correlation coefficient for all conditions
corr_result = cor(expression_with_bem_mat %>% select(-index, -ensembl_gene_id, -gene_id, -gene_symbol, -entrez_id), method = "spearman", use = "complete.obs")
corr_result_df = as.tibble(corr_result)

# plot heatmap 
png('~/Projects/yonsei_organoid/analysis/corr_coeff_with_brainspan.png')
pheatmap(corr_result)
dev.off() 



corr_result_df$index = row.names(corr_result_df)

corr_result_df %>% gather(c(1:524, 'mat', 'bem', 'tissue', 'nsc'), key='column', value='corr_coef') %>% filter(corr_coef <1, column=='mat')
View(corr_result_df %>% gather(c(1:524, 'mat', 'bem'), key='column', value='corr_coef') %>% filter(corr_coef <1. &  column=='bem') %>% arrange(desc(corr_coef)))

#plot correlation coefficients
png('~/Projects/yonsei_organoid/analysis/mat_spearman_coeff_hist.png')
mat_coefficients = (corr_result_df %>% gather(c(1:524, 'mat', 'bem', 'tissue', 'nsc'), key='column', value='corr_coef') %>% filter(corr_coef <1. &  column=='mat', index %in% 1:524) %>% arrange(desc(corr_coef)))
mat_mean_coefficient = as.double(mat_coefficients %>% summarise(mean=mean(corr_coef)))
ggplot(mat_coefficients, aes(corr_coef)) + geom_histogram() + xlim(0.5,1) + ylim(0, 400) + ggtitle('MAT spearman coefficients with Brainspan RNA-Seq Data') + geom_vline(xintercept = mat_mean_coefficient) + annotate('text', x=mat_mean_coefficient, y= 300, label=str_c('mean coefficient =', round(mat_mean_coefficient, digits = 4)))
dev.off() 

png('~/Projects/yonsei_organoid/analysis/bem_spearman_coeff_hist.png')
bem_coefficients = (corr_result_df %>% gather(c(1:524, 'mat', 'bem', 'tissue', 'nsc'), key='column', value='corr_coef') %>% filter(corr_coef <1. &  column=='bem', index %in% 1:524) %>% arrange(desc(corr_coef)))
bem_mean_coefficient = as.double(bem_coefficients %>% summarise(mean=mean(corr_coef)))
ggplot(bem_coefficients, aes(corr_coef)) + geom_histogram() + xlim(0.5,1) + ylim(0, 400) + ggtitle('BEM spearman coefficients with Brainspan RNA-Seq Data')+ geom_vline(xintercept = bem_mean_coefficient) + annotate('text', x=bem_mean_coefficient, y= 300, label=str_c('mean coefficient =', round(bem_mean_coefficient, digits = 4)))
dev.off() 

png('~/Projects/yonsei_organoid/analysis/tissue_spearman_coeff_hist.png')
tissue_coefficients = (corr_result_df %>% gather(c(1:524, 'mat', 'bem', 'tissue', 'nsc'), key='column', value='corr_coef') %>% filter(corr_coef <1. &  column=='tissue', index %in% 1:524) %>% arrange(desc(corr_coef))) 
tissue_mean_coefficient = as.double(tissue_coefficients %>% summarise(mean=mean(corr_coef)))
ggplot(tissue_coefficients, aes(corr_coef)) + geom_histogram() + xlim(0.5,1)+ ylim(0, 400)  + ggtitle('Tissue spearman coefficients with Brainspan RNA-Seq Data')+ geom_vline(xintercept = tissue_mean_coefficient) + annotate('text', x=tissue_mean_coefficient, y= 300, label=str_c('mean coefficient =', round(tissue_mean_coefficient, digits = 4)))
dev.off()

png('~/Projects/yonsei_organoid/analysis/nsc_spearman_coeff_hist.png')
nsc_coefficients = (corr_result_df %>% gather(c(1:524, 'mat', 'bem', 'tissue', 'nsc'), key='column', value='corr_coef') %>% filter(corr_coef <1. &  column=='nsc', index %in% 1:524) %>% arrange(desc(corr_coef)))
nsc_mean_coefficient = as.double(nsc_coefficients %>% summarise(mean=mean(corr_coef)))
ggplot(nsc_coefficients, aes(corr_coef)) + geom_histogram() + xlim(0.5,1) + ylim(0, 400) + ggtitle('NSC spearman coefficients with Brainspan RNA-Seq Data')+ geom_vline(xintercept = nsc_mean_coefficient) + annotate('text', x=nsc_mean_coefficient, y= 300, label=str_c('mean coefficient =', round(nsc_mean_coefficient, digits = 4)))
dev.off() 
