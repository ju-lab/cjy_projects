# R script to generate Immunoglobulin Gene Locuss with V, (D), J, and C regions 
# and mark locations of structural variations in the same plot 
# 2018.07.20 CJY

# library(grid)
library(tidyverse)
library(vcfR)
setwd("~/Documents/cjyoon/Projects/myeloma/analysis/ig")
ig_genes = as.tibble(read.table('immunoglobulin_grch37.bed', col.names =  c("chromosome", "start", "end", "gene",  "strand", 'type')))
ig_genes$strand <- ifelse(ig_genes$strand == "+", 'forward', 'reverse')
ig_genes

igh = ig_genes %>% dplyr::filter(chromosome=='chr14')
igl = ig_genes %>% dplyr::filter(chromosome =='chr22')
igk = ig_genes %>% dplyr::filter(chromosome=='chr2')

# Function to draw each genes
draw_gene <- function(gene_bed){
  min_pos = min(gene_bed$start)
  max_pos= max(gene_bed$end)
  max_dist = max_pos -min_pos
  plot_min = min_pos - 0.1 * max_dist
  plot_max = max_pos + 0.1 * max_dist
  plot(0, type = "n", xlim=c(plot_min, plot_max))
  title(unique(gene_bed$chromosome))
  gene_bed$median = 0.5 * (gene_bed$start + gene_bed$end)
  
  gene_bed$color = gene_bed$type %>% str_replace(".*-", "") %>% factor(levels=c("V", "D", "J", "C")) %>% as.numeric %>% {rainbow(4)[.]}
  with(gene_bed, rect(start, -0.4, end, -0.1, col = color, lty=0.01))
  with(gene_bed, points(median, rep(0, nrow(gene_bed)), col = color, pch=20 ))
}
draw_gene2 <- function(gene_bed, level = -1, max_level = -1){
  min_pos = min(gene_bed$start)
  max_pos= max(gene_bed$end)
  max_dist = max_pos -min_pos
  plot_min = min_pos - 0.1 * max_dist
  plot_max = max_pos + 0.1 * max_dist
  plot_r = plot_max - plot_min
  plot(0, type = "n", xlim=c(0, 1))
  gene_bed$start.f <- (gene_bed$start - plot_min) / plot_r
  gene_bed$end.f <- (gene_bed$end - plot_min) / plot_r
  title(unique(gene_bed$chromosome))
  gene_bed$median = 0.5 * (gene_bed$start.f + gene_bed$end.f)
  
  gene_bed$color = gene_bed$type %>% str_replace(".*-", "") %>% factor(levels=c("V", "D", "J", "C")) %>% as.numeric %>% {rainbow(4)[.]}
  with(gene_bed, rect(start.f, -0.4, end.f, -0.1, col = color, lty=0.01))
  with(gene_bed, points(median, rep(0, nrow(gene_bed)), col = color, pch=20 ))
  pretty(c(plot_min, plot_max)) %>% print() %>% .[2:(length(.)-1)] %>%{text(x = (. - plot_min)/plot_r, y = -0.5, .)}
  
}
par(mfrow=c(2,1))
draw_gene2(igl)
draw_gene(igl)

# Function to draw each genes into separate graphs
draw_all_genes <- function(gene_bed){
  nchrom =unique(gene_bed$chromosome)
  par(mfrow=c(length(nchrom), 1))
  for(chrom in nchrom){
    gene_bed %>% dplyr::filter(chromosome==chrom) %>% draw_gene
  }
  legend("bottomright", inset=.02, title="IG segments",
         c("V", "D", "J", "C"), fill=rainbow(4), horiz=TRUE, cex=0.8)
}

pdf('ig_genes.pdf')
draw_all_genes(ig_genes)
dev.off()

# Function to draw a curve between two points (x1, y1) and (x2, y2)
curve <- function(x1, y1, x2, y2, color, scale){
  X<- seq(x1,x2,length.out=100)
  Y<- seq(y1,y2,length.out=100)
  C<-dbeta(1:100/100,2,2)
  # scale = max(abs(y2-y1), 2) * 0.1
  lines(X,Y+C*scale, col = color)
}

transform<-function(x, y, scale, theta, h, v){
  x_t = scale* ( x*cos(theta) - y * sin(theta)) + h
  y_t = (x*sin(theta) + y * cos(theta)) + v
  
  return(list(x_t, y_t))
}


sine_curve <- function(x1, y1, x2, y2, extra_scale=1, color='black'){
  # Y=sin(X)
  X<- seq(0, pi, length.out=100)
  Y<- extra_scale * sin(X)
  
  # length scale
  scale = sqrt((y2-y1)**2 + (x2-x1)**2) / (pi)
  # Rotation angle
  theta = atan((y2-y1)/(x2-x1))
  
  # transformed vectors
  transformed_sin = transform(X, Y, scale, theta, min(x1,x2), min(y1,y2))
  lines(transformed_sin[1][[1]],transformed_sin[2][[1]], col = color)
  return(transformed_sin)
}



draw_gene_locus <- function(gene_bed, locus){
  min_pos = min(gene_bed$start)
  max_pos= max(gene_bed$end)
  max_dist = max_pos -min_pos
  plot_min = min_pos - 0.1 * max_dist
  plot_max = max_pos + 0.1 * max_dist
  plot(0, type = "n", xlim=c(plot_min, plot_max))
  title(unique(gene_bed$chromosome))
  gene_bed$median = 0.5 * (gene_bed$start + gene_bed$end)
  gene_bed$color = gene_bed$type %>% str_replace(".*-", "") %>% factor(levels=c("V", "D", "J", "C")) %>% as.numeric %>% {rainbow(4)[.]}
  with(gene_bed, rect(start, -0.4, end, -0.1, col = color, lty=0.01))
  with(gene_bed, points(median, rep(0, nrow(gene_bed)), col = color, pch=20 ))
  with(del_locus, points(start1, 0.5, start2, 0.5))
  
}
 

##################
# Parse VCF file
complement_base = c('A', 'T', 'G', 'C')
names(complement_base) = c('T', 'A', 'C', 'G')

# convert 12 type base changes into 6 type
basechange <- function(ref, alt){
  ifelse(ref %in% c("A", "G"), paste0(complement_base[ref], '>', complement_base[alt]),paste0(ref, '>', alt) )
}



# including absCN and SNVs like Yilong Plot
gene_with_sv_locus <- function(gene_bed, locus, absCN, snv_vcf, legend=T){
  min_pos = min(gene_bed$start)
  max_pos= max(gene_bed$end)
  max_dist = max_pos -min_pos
  plot_min = min_pos - 0.1 * max_dist
  plot_max = max_pos + 0.1 * max_dist
  plot_r = plot_max - plot_min
  
  gene_chromosome = unique(gene_bed$chromosome)
  
  
  plot(0, type = "n", bty="n", ylab=gene_chromosome, xlab='', yaxt='n', axes=FALSE,  xlim=c(0, 1), ylim=c(-1, 2))
  if(gene_chromosome=='chr14'){
    title('IGH')
  }else if(gene_chromosome=='chr22'){
    title('IGL')
  }else if(gene_chromosome=='chr2'){
    title('IGK')
  }else{
    title('Unknown')
  }
  gene_bed$start.f <- (gene_bed$start - plot_min) / plot_r
  gene_bed$end.f <- (gene_bed$end - plot_min) / plot_r
  
  gene_chromosome_number = gsub('chr', '', gene_chromosome)
  gene_bed$median = 0.5 * (gene_bed$start.f + gene_bed$end.f)
  
  gene_bed$color = gene_bed$type %>% str_replace(".*-", "") %>% factor(levels=c("V", "D", "J", "C")) %>% as.numeric %>% {rainbow(4)[.]}
  with(gene_bed, rect(start.f, -0.4, end.f, -0.1, col = color, lty=0.01))
  with(gene_bed, points(median, rep(0, nrow(gene_bed)), col = color, pch=20 ))
  
  
  # prepare locus dataframe, only if locus_df is not an empty data file
  if(dim(locus)[1]!=0){
    # separate position into separate columns with chromosome/start/end
    locus_df = locus %>% separate(bp1, into=c('bp1_chr', 'bp1_start', 'bp1_end'), sep='[:-]') %>% separate(bp2, into=c('bp2_chr', 'bp2_start', 'bp2_end'), sep='[:-]') %>%  mutate(keep_forplot = ((bp1_chr) == gene_chromosome_number | bp2_chr == gene_chromosome_number))
    
    # to avoid positions in other locus from being plotted here. 
    locus_df = locus_df %>% filter(keep_forplot==TRUE)
    locus_df$bp1_start = as.integer(locus_df$bp1_start)
    locus_df$bp2_start = as.integer(locus_df$bp2_start)
    
    locus_df$start.f <- (locus_df$bp1_start - plot_min)/plot_r
    locus_df$end.f <- (locus_df$bp2_start - plot_min)/plot_r
    
    bnd_locus_df = locus_df %>% dplyr::filter(svtype=='BND')
    del_locus_df = locus_df %>% dplyr::filter(svtype=='DEL')
    inv_locus_df = locus_df %>% dplyr::filter(svtype=='INV')
    dup_locus_df = locus_df %>% dplyr::filter(svtype=='DUP')
    
    # Plot Deletion  
    del_pos = 0.3
    segments(0, del_pos, 1, del_pos, col='gray')
    text(0, del_pos+0.15, "DEL", col = "gray", cex=0.8)
    if(dim(del_locus_df)[1]>0){
      for(i in 1:dim(del_locus_df)[1]){
        sine_curve(del_locus_df[i, 'start.f'][[1]], del_pos, del_locus_df[i, 'end.f'][[1]], del_pos, extra_scale=0.2, color='black')
      }
    }
    dup_pos = 0.3
    text(0, dup_pos-0.15, "TD", col = "gray", cex=0.8)
    
    if(dim(dup_locus_df)[1] >0){
      for(i in 1:dim(dup_locus_df)[1]){
        sine_curve(dup_locus_df[i, 'start.f'][[1]], dup_pos, dup_locus_df[i, 'end.f'][[1]], dup_pos, extra_scale=-0.2, color='black')
      }
    }
    
    
    inv_pos = 0.8
    segments(0, inv_pos, 1, inv_pos, col='gray')
    text(0, inv_pos+0.15, "TT", col = "gray", cex=0.8)
    text(0, inv_pos-0.15, "HH", col = "gray", cex=0.8)
    
    if(dim(inv_locus_df)[1]>0){
      for(i in 1:dim(inv_locus_df)[1]){
        if(inv_locus_df[i, 'connection'] =='3to3'){ # tail to tail 
          sine_curve(inv_locus_df[i, 'start.f'][[1]], inv_pos, inv_locus_df[i, 'end.f'][[1]], inv_pos, extra_scale=0.2, color='black')
        }else if(inv_locus_df[i, 'connection'] == '5to5'){
          sine_curve(inv_locus_df[i, 'start.f'][[1]], inv_pos, inv_locus_df[i, 'end.f'][[1]], inv_pos, extra_scale=-0.2, color='black')
          
        }
        
      }
    }
    tra_pos = 1.3
    if(dim(bnd_locus_df)[1] > 0){
      for(i in 1:dim(bnd_locus_df)[1]){
        if(bnd_locus_df[i, 'connection'] == '5to5'){
          label_xpos = bnd_locus_df[i, 'start.f'][[1]]+0.05
          label_ypos = tra_pos + 0.5
        }else if(bnd_locus_df[i, 'connection'] == '5to3'){
          label_xpos = bnd_locus_df[i, 'start.f'][[1]] - 0.05
          label_ypos = tra_pos + 0.5
        }else if(bnd_locus_df[i, 'connection'] == '3to5'){
          label_xpos = bnd_locus_df[i, 'start.f'][[1]]+0.05
          label_ypos = tra_pos + 0.5
        }else if(bnd_locus_df[i, 'connection'] == '3to3'){
          label_xpos = bnd_locus_df[i, 'start.f'][[1]] - 0.05
          label_ypos = tra_pos + 0.5          
        }else{
          print(0)
        }
        curve(bnd_locus_df[i, 'start.f'][[1]], tra_pos, label_xpos, label_ypos, scale=0.3, color='black')
        breakpoint_chromosomes = c(bnd_locus_df[i,'bp1_chr'][[1]], bnd_locus_df[i,'bp2_chr'])
        translocated_chromosome = breakpoint_chromosomes[breakpoint_chromosomes != gene_chromosome_number]
        text(label_xpos+0.02, label_ypos, paste('chr', translocated_chromosome, sep=''))
      }
    }
  }
  # 
  # with(del_locus_df, grid.curve(start.f, rep(0.3, length(del_locus_df$end.f)), end.f, rep(0.3, length(del_locus_df$end.f)), curvature=1))
  # with(inv_locus_df, grid.curve(start.f, rep(0.4, length(inv_locus_df$end.f)), end.f, rep(0.4, length(inv_locus_df$end.f))))
  # 
  
  # Plot intermutation distances in this region
  vcf <- read.vcfR(snv_vcf)
  vcf_tibble <- as.tibble(getFIX(vcf))
  vcf_tibble$POS = as.integer(vcf_tibble$POS)
  vcf_logdist = vcf_tibble %>% mutate(basechange = basechange(REF, ALT)) %>% mutate(distance = POS- lag(POS)) %>% mutate(log_dist = log10(distance)) %>% filter(CHROM==str_remove(gene_chromosome, 'chr') & between(POS, min_pos, max_pos))
  vcf_logdist$color = vcf_logdist$basechange %>% factor(levels=c("T>A", "T>G", "T>C", "C>A", "C>G", "C>T")) %>% as.numeric %>% {rainbow(6)[.]}
  vcf_logdist$POS.f <- (vcf_logdist$POS - plot_min) / plot_r
  plot(vcf_logdist$POS.f, vcf_logdist$log_dist, col=vcf_logdist$color, xlim=c(0, 1), ylim=c(0, 6),  bty="n", ylab='log10(mutation distance)', xlab='', xaxt='n')
  
  # Plot absCN in this region
  absCN_df = read_delim(absCN, skip=4, delim='\t', col_names=c('chrom', 'pos', 'dummy1', 'dummy2', 'abscn'))
  absCN_locus_df = absCN_df %>% filter(chrom==str_remove(gene_chromosome, 'chr') & between(pos, min_pos, max_pos))
  absCN_locus_df$pos.f = (absCN_locus_df$pos - plot_min) / plot_r
  plot(absCN_locus_df$pos.f, absCN_locus_df$abscn, xlim=c(0, 1), ylim=c(0, 6),  bty="n", ylab='absCN', xlab='', xaxt='n')
  

  pretty(c(plot_min, plot_max)) %>% print() %>% .[2:(length(.)-1)] %>%{text(x = (. - plot_min)/plot_r, y = -0.5, .)}
  if(legend==TRUE){
    legend("topright", inset=.02, title="IG segments",
           c("V", "D", "J", "C"), fill=rainbow(4), horiz=TRUE, cex=0.8)
  }
}


# Actual part of the code that reads in SV calls and plots SVs in the IG regions

files <-  Sys.glob(file.path('../sv_intersect/', "*_sv_filtered.intersectFilter.vcf.svanno.txt"))
for(file in files){
  locus = read_delim(file, col_names=c('bp1', 'bp2', 'svtype', 'connection', 'bp1annot', 'bp2_annot'), delim='\t')
  samplename = gsub('../sv_intersect//', '', file)
  samplename = gsub('_sv_filtered.intersectFilter.vcf.svanno.txt', '', samplename)
  print(samplename)
  pdf(paste(samplename, '_IG_SVs_v2.pdf', sep=''), width=7, height=20)
  
  snv_vcf = paste0('~/Documents/cjyoon/Projects/myeloma/analysis/snv_filter_v2/', samplename, '_mns.pon.filtered.vcf')
  absCN = paste0('~/Documents/cjyoon/Projects/myeloma/analysis/abscn/', samplename, '.100kbcov.absCN')
  
  par(mfrow=c(9,1), xpd=NA, oma = c(4, 4, 2, 0), mar=c(4,1,1,1))
  gene_with_sv_locus(igh, locus, absCN, snv_vcf, legend=FALSE)
  gene_with_sv_locus(igk, locus, absCN, snv_vcf, legend=FALSE)
  gene_with_sv_locus(igl, locus,absCN, snv_vcf, legend=TRUE)
  mtext(samplename, outer = TRUE, cex = 1.5)
  dev.off()
  
}




snv_vcf = '~/Documents/cjyoon/Projects/myeloma/analysis/snv_filter_v2/3346_mns.pon.filtered.vcf'
absCN =  '~/Documents/cjyoon/Projects/myeloma/analysis/abscn/3346.100kbcov.absCN'

gene_with_sv_locus(igh, locus, absCN = absCN , snv_vcf = snv_vcf, legend=FALSE)
dev.off()
ggsave("~/Documents/dist.png")
