# This script plots large deletion duplication identified in PCAWG mitochondira
# Figure 5G of the paper
# April 5 2018, Chris Yoon (cjyoon@kaist.ac.kr)

library(tidyverse)
library(ggplot2)
library(scater)
library(gggenes)
library(ggpubr)
library(RcppRoll)
library(cowplot)
# Import per-base depth data 
sample_1_cov_data = '~/Downloads/12618023223ce42839d65965b478c4c6.bam.depth.MT'
sample_1_cov_data_df = read_delim(sample_1_cov_data , delim = '\t', col_names = c('chromosome', 'pos', 'depth'))

sample_2_cov_data = '~/Downloads/ac570307b3acc0cdc0a250cd1f13ada2.bam.depth.MT'
sample_2_cov_data_df = read_delim(sample_2_cov_data, delim = '\t', col_names = c('chromosome', 'pos', 'depth'))

sample_3_cov_data = '~/Downloads/49a5f51897e72742aa52eaf1134205b0.bam.depth.MT'
sample_3_cov_data_df = read_delim(sample_3_cov_data, delim = '\t', col_names = c('chromosome', 'pos', 'depth'))


# Define CNV regions based on observation 
sample_1_cnv <- function(position){ifelse(between(position, 13464,14891), 'del', 'normal')}
sample_2_cnv <- function(position){ifelse(between(position, 755,3557), 'del', 'normal')}
sample_3_cnv <- function(position){ifelse(between(position, 12014,16270), 'dup', 'normal')}


# calculate average depth for normal and cnv regions for all samples
s1_avg_cnv=sample_1_cov_data_df %>% dplyr::mutate(is_cnv = sample_1_cnv(pos)) %>% group_by(is_cnv) %>% summarise(mean=mean(depth))
s2_avg_cnv=sample_2_cov_data_df %>% dplyr::mutate(is_cnv = sample_2_cnv(pos)) %>% group_by(is_cnv) %>% summarise(mean=mean(depth))
s3_avg_cnv=sample_3_cov_data_df %>% dplyr::mutate(is_cnv = sample_3_cnv(pos)) %>% group_by(is_cnv) %>% summarise(mean=mean(depth))

# Calculate percentage of deleted/duplicated mitochondria
s1_percent = 1- s1_avg_cnv[s1_avg_cnv$is_cnv=='del', 'mean']/s1_avg_cnv[s1_avg_cnv$is_cnv=='normal', 'mean']
s2_percent = 1- s2_avg_cnv[s2_avg_cnv$is_cnv=='del', 'mean']/s2_avg_cnv[s2_avg_cnv$is_cnv=='normal', 'mean']
s3_percent = s3_avg_cnv[s3_avg_cnv$is_cnv=='dup', 'mean']/s3_avg_cnv[s3_avg_cnv$is_cnv=='normal', 'mean'] -1 


# Mitochondria gene annotation

mito_genes = as.tibble(read.table('/Users/cyoon/Documents/julab/projects/mito/mito_genes.bed', col.names =  c("molecule", "start", "end", "gene", "dummy",  "strand", 'type')))
mito_genes$direction <- ifelse(mito_genes$strand == "+", 1, -1)

# Plot
pdf('~/Downloads/pcawg_mito_deldup_v3.pdf')
margin(t = 0, r = 0, b = 0, l = 0, unit = "cm")


s1 <- sample_1_cov_data_df %>% dplyr::mutate(rolling_mean = roll_mean(x=depth,n=100, align='center', na.rm = T, fill=NA)) %>% 
  ggplot(aes(x=pos, y=rolling_mean)) + geom_point() + ggtitle(str_c('SA501399 LIHC'))+
  scale_x_continuous(limits = c(0,16569), expand = c(0, 0),  breaks = seq(0, 16000, by = 4000)) + scale_y_continuous(limits = c(0,max(sample_1_cov_data_df$depth))) + 
  geom_segment(aes(x=13464, y=28774, xend=14891, yend=28774), color='green') +
  geom_segment(aes(x=0, y=43900, xend=13463, yend=43900), color='grey', linetype=2) + 
  geom_segment(aes(x=14892, y=43900, xend=16569, yend=43900), color='grey', linetype=2)+
  xlab('') + ylab('Depth') + theme_pubr() + theme(axis.title =element_text(size=0), aspect.ratio = 0.15, plot.margin = margin(-3, 0, -3, 0, "cm"), axis.title.x=element_blank(), panel.grid = element_blank(), panel.border = element_blank()) +  annotate("text", Inf, 10000, label = str_c('Deletion in ', round(s1_percent, digits = 2)*100, '%'), hjust='top', vjust = 'right', size=3, colour='grey')
ggsave('~/Downloads/SA501399.mitodel.pdf',  width = 10, height = 2, dpi = 150)

s2 <-sample_2_cov_data_df %>% dplyr::mutate(rolling_mean = roll_mean(x=depth,n=100, align='center', na.rm = T, fill=NA)) %>% 
  ggplot(aes(x=pos, y=rolling_mean)) + geom_point() + ggtitle(str_c('SP76017 PACA'))+
  scale_x_continuous(limits = c(0,16569), expand = c(0, 0),  breaks = seq(0, 16000, by = 4000))+ scale_y_continuous(limits = c(0,max(sample_2_cov_data_df$depth))) +
  geom_segment(aes(x=755, y=3327, xend=3557, yend=3327), color='green') +
  geom_segment(aes(x=0, y=8891, xend=754, yend=8891), color='grey', linetype=2) + 
  geom_segment(aes(x=3558, y=8891, xend=16569, yend=8891), color='grey', linetype=2)+
xlab('') + ylab('Depth')+ theme_pubr() + theme(aspect.ratio = 0.15, plot.margin = margin(-3, 0, -3, 0, "cm"), axis.title.x=element_blank(), panel.grid = element_blank(), panel.border = element_blank()) + annotate("text", 4000, 3000, label=str_c('Deletion in ', round(s2_percent, digits = 2)*100, '%'), hjust='left', vjust = 'right', size=3, colour='grey') 
ggsave('~/Downloads/SP76017.mitodel.pdf',  width = 10, height = 2, dpi = 150)

s3<-sample_3_cov_data_df %>% dplyr::mutate(rolling_mean = roll_mean(x=depth,n=100, align='center', na.rm = T, fill=NA)) %>% 
  ggplot(aes(x=pos, y=rolling_mean)) + geom_point() + ggtitle('SP127680 MELA')+
  scale_x_continuous(limits = c(0,16569), expand = c(0, 0), breaks = seq(0, 16000, by = 4000)) + scale_y_continuous(limits = c(0,max(sample_3_cov_data_df$depth))) + 
  geom_segment(aes(x=12014, y=34197, xend=16270, yend=34197), color='red') + 
  geom_segment(aes(x=0, y=15230, xend=12013, yend=15230), color='grey', linetype=2) + 
  geom_segment(aes(x=16271, y=15230, xend=16569, yend=15230), color='grey', linetype=2) +
  xlab('MT position') + ylab('Depth')+ theme_pubr() + theme(aspect.ratio = 0.15, plot.margin = margin(-3, 0, -3, 0, "cm"), axis.title.x=element_blank(), panel.grid = element_blank(), panel.border = element_blank())+ annotate("text", 13000, 10000, label=str_c('Duplication'), hjust='left', vjust = 'right', size=3, colour='grey')+  
  theme(plot.margin = unit(c(0,0,0,0), "cm"))
ggsave('~/Downloads/SP127680.mitodup.pdf', width = 10, height = 2, dpi = 150)

genes <- ggplot2::ggplot(mito_genes, ggplot2::aes(xmin = start, xmax = end, y =
                                                    molecule, fill = type, label = gene, forward = direction)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")) +
  geom_gene_label(align = "left") +
  ggplot2::facet_wrap(~ molecule, scales = "free", ncol = 1) +
  ggplot2::scale_fill_brewer(palette = "Set3") +
  theme_genes() + guides(fill=FALSE)+ theme(legend.position="bottom") + ylab('') + xlab('MT position')+ scale_x_continuous(limits = c(0,16569), expand = c(0, 0), breaks = seq(0, 16000, by = 4000))#, axis.title.y=element_blank())
# grid.draw(rbind(ggplotGrob(s1), ggplotGrob(s2), ggplotGrob(s3), ggplotGrob(genes), size = "last"))
# ggarrange(s1, s2, s3, genes, heights = c(10, 20),
          # ncol = 1, nrow = 4, align = "h")
ggsave('~/Downloads/mito_genes.pdf', width = 10, height = 1, dpi = 150)
# 
plot_grid(s1, s2, s3, genes,  align = 'v', axis = 'l', ncol=1, labels=c('G', '', '', ''), rel_heights = c(0.5, 0.5,0.5,0.2))
dev.off()




