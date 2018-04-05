library(tidyverse)
library(ggplot2)
library(scater)
setwd('~/Documents/julab/projects/pcawg/Archive/')

tumor_depth_concat = '~/Documents/julab/projects/pcawg/Archive/pcawg_tumourMT.100bpbin.cov'
normal_depth_concat = '~/Documents/julab/projects/pcawg/Archive/pcawg_normalMT.100bpbin.cov' 


# This function is to get normalized count per sample. 
normalized_depth <- function(depth_concat){
  data = read.table(depth_concat, col.names=c('bin', 'mean_cov', 'bam'))
  data = as.tibble(data)
  mean_depth = data %>% dplyr::group_by(bam) %>% dplyr::summarise(mean = mean(mean_cov), std = sd(mean_cov))
  normalized = data %>% dplyr::left_join(mean_depth, by='bam') 
  normalized %>% dplyr::mutate(normalized = (mean_cov/mean))
  
}


# Use the function above to get normlized count within each sample (divided by average depth)
tumor_normalized = normalized_depth(tumor_depth_concat)
normal_normalized = normalized_depth(normal_depth_concat)

tumor_normalized$type = 't'
normal_normalized$type = 'n'

# Compare each nomral against other normals
# Create a funciton to create a plot for each normal sample against panel of normal (excluding itself)

plot_eachnormalbam <- function(normal_normalized, normal_sample){
  print(normal_sample)
  normal_normalized_minus_itself = normal_normalized %>% dplyr::filter(bam != normal_sample)
  normal_normalized_itself =  normal_normalized %>% dplyr::filter(bam == normal_sample)
  normal_normalized_by_pos_minus_itself = normal_normalized_minus_itself %>% dplyr::group_by(bin) %>% dplyr::summarise(pos_normalized_mean=mean(normalized), pos_normalized_std = sd(normalized))
  normal_normalized_pon_minus_itself = normal_normalized_itself %>% dplyr::left_join(normal_normalized_by_pos_minus_itself, by='bin') %>% dplyr::mutate(normalized_pon = normalized/pos_normalized_mean) %>% dplyr::mutate(zscore = (normalized-pos_normalized_mean)/pos_normalized_std)
  ggplot(normal_normalized_pon_minus_itself, aes(x=bin, y=zscore, color=bam)) + geom_line() + scale_y_continuous(limits = c(-10,10)) + theme(legend.position="none") + ggtitle(normal_sample)
  
}

normal_bamlist = unique(normal_normalized$bam)

for(i in 1:as.integer(length(normal_bamlist)/8)){
  bamnumber =(8*i-7):(8*i)
  p1 = plot_eachnormalbam(normal_normalized, normal_bamlist[bamnumber[1]])
  p2 = plot_eachnormalbam(normal_normalized, normal_bamlist[bamnumber[2]])
  p3 = plot_eachnormalbam(normal_normalized, normal_bamlist[bamnumber[3]])
  p4 = plot_eachnormalbam(normal_normalized, normal_bamlist[bamnumber[4]])
  p5 = plot_eachnormalbam(normal_normalized, normal_bamlist[bamnumber[5]])
  p6 = plot_eachnormalbam(normal_normalized, normal_bamlist[bamnumber[6]])
  p7 = plot_eachnormalbam(normal_normalized, normal_bamlist[bamnumber[7]])
  p8 = plot_eachnormalbam(normal_normalized, normal_bamlist[bamnumber[8]])
  pdf(str_c('~/Documents/julab/projects/pcawg/Archive/normal_plots/Normal_batchnumber', i, '.pdf', sep=''),  width=12, height=14)
  scater::multiplot(p1, p2, p3, p4, p5, p6, p7, p8, cols=2)
  dev.off()
}





# Create combined tibble to compare tumor against panel of normal
combined = bind_rows(tumor_normalized, normal_normalized)
ggplot(combined, aes(x=bin, y=normalized, colour=bam, shape=as.factor(type))) + geom_point() + theme(legend.position="none")



normal_normalized_by_pos = normal_normalized %>% group_by(bin) %>% summarise(pos_normalized_mean=mean(normalized), pos_normalized_std = sd(normalized))
tumor_normalized_pon = tumor_normalized %>% left_join(normal_normalized_by_pos, by='bin') %>% mutate(normalized_pon = normalized/pos_normalized_mean) %>% mutate(zscore = (normalized-pos_normalized_mean)/pos_normalized_std)
ggplot(tumor_normalized_pon, aes(x=bin, y=zscore, color=bam)) + geom_point() + theme(legend.position="none")


plot_eachbam <- function(tumor_normalized_pon, bamfile){
  tumor_normalized_pon %>% filter(bam==bamfile) %>% ggplot(aes(x=bin, y=zscore)) + geom_line() + scale_y_continuous(limits = c(-10,10)) + ggtitle(bamfile)
}

bamlist = unique(tumor_normalized_pon$bam)

for (i in 1:as.integer(length(bamlist)/8)){
  bamnumber =(8*i-7):(8*i)
  print(i)
  p1 = plot_eachbam(tumor_normalized_pon, bamlist[bamnumber[1]])
  p2 = plot_eachbam(tumor_normalized_pon, bamlist[bamnumber[2]])
  p3 = plot_eachbam(tumor_normalized_pon, bamlist[bamnumber[3]])
  p4 = plot_eachbam(tumor_normalized_pon, bamlist[bamnumber[4]])
  p5 = plot_eachbam(tumor_normalized_pon, bamlist[bamnumber[5]])
  p6 = plot_eachbam(tumor_normalized_pon, bamlist[bamnumber[6]])
  p7 = plot_eachbam(tumor_normalized_pon, bamlist[bamnumber[7]])
  p8 = plot_eachbam(tumor_normalized_pon, bamlist[bamnumber[8]])
  pdf(str_c('~/Documents/julab/projects/pcawg/Archive/plots/batchnumber', i, '.pdf', sep=''),  width=12, height=14)
  scater::multiplot(p1, p2, p3, p4, p5, p6, p7, p8, cols=2)
  dev.off()
    }
plot_eachbam(tumor_normalized_pon, bamlist[20])

# Look at those with weird coverage profile
weird_list = c('0095dee3b343c04564a762b9ff2dd316.bam', '12618023223ce42839d65965b478c4c6.bam', '25159a273d764e2c4dfad79cc655217e.bam', '28fd6d07302e17e9a60a2dff0a0ed9da.bam', '47e0c917aecce2affd8bf4c696a8eb5e.bam', '3bcb0039e46a4ca217d29b0c77a1ef7e.bam', '3c17cc820c10e38748d50dec11367cfb.bam', '49a5f51897e72742aa52eaf1134205b0.bam', '663f47b32ec2485f698c5cd2110c4ca1.bam', '7167f0e3e4ae49d781d7204718e164e2.bam', '7fd88780cf59cfb445801dd869e84a40.bam', '9aeba91791085198536ec317e18eecd9.bam', 'a539dd7d8277fdc875427473edc12e7c.bam', 'ac570307b3acc0cdc0a250cd1f13ada2.bam', 'd1dc0e2c2afe463f2cbf21ff03dbe746.bam', 'e5831d43ef17ecd0273a8ba5d3647ebc.bam')

for (i in 1:as.integer(length(weird_list)/8)){
  bamnumber =(8*i-7):(8*i)
  print(i)
  p1 = plot_eachbam(tumor_normalized_pon, weird_list[bamnumber[1]])
  p2 = plot_eachbam(tumor_normalized_pon, weird_list[bamnumber[2]])
  p3 = plot_eachbam(tumor_normalized_pon, weird_list[bamnumber[3]])
  p4 = plot_eachbam(tumor_normalized_pon, weird_list[bamnumber[4]])
  p5 = plot_eachbam(tumor_normalized_pon, weird_list[bamnumber[5]])
  p6 = plot_eachbam(tumor_normalized_pon, weird_list[bamnumber[6]])
  p7 = plot_eachbam(tumor_normalized_pon, weird_list[bamnumber[7]])
  p8 = plot_eachbam(tumor_normalized_pon, weird_list[bamnumber[8]])
  pdf(str_c('~/Documents/julab/projects/pcawg/Archive/plots/weird_bams_', i, '.pdf', sep=''),  width=12, height=14)
  scater::multiplot(p1, p2, p3, p4, p5, p6, p7, p8, cols=2)
  dev.off()
}

