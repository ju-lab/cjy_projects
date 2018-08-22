library(tidyverse)
library(stringr)
library(biomaRt)
library(refGenome)

ig_annotation_path = '/Users/cyoon/Documents/cjyoon/Projects/myeloma/analysis/ig/immunoglobulin_grch37.bed'
ig_annotation_df = read_delim(ig_annotation_path, delim='\t', col_names = c("chromosome", "start", "end", "symbol", "strand", "classification"))
ig_annotation_df = ig_annotation_df %>% separate(classification, into=c('chain', 'VDJC'))

# # get Ensembl to HUGO gene symbol table
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mapping <- as.tibble(getBM(attributes = c("ensembl_gene_id", 'hgnc_symbol', 'chromosome_name', 'start_position', 'end_position'), mart = ensembl))

# For each RSEM output draw log2 (TPM + 1) value of immunoglobulin gene segments. 

plot_ig_expressions <- function(sampleID){
  star_rsem_path = paste0('~/Documents/cjyoon/scripts/star/rsem/', sampleID, '.genes.results')
  new_tuxedo_path = paste0('/Users/cyoon/Documents/cjyoon/scripts/new_tuxedo/stringtie/quant/', sampleID, '/', sampleID, '.transcript.gtf')
  
  star_rsem_df = read_delim(star_rsem_path, delim='\t')
  
  
  # star_rsem_df
  # 
  # # 
  ig_expressions = star_rsem_df %>% left_join(mapping, by=c('gene_id'='ensembl_gene_id')) %>% filter(chromosome_name %in% c('14', '2', '22')) %>% filter(str_detect(hgnc_symbol, '^IGH|^IGL|^IGK')) %>% filter(!str_detect(hgnc_symbol, 'OR')) %>% arrange(chromosome_name, start_position)
  ig_expressions = ig_expressions %>% left_join(ig_annotation_df, by=c('hgnc_symbol'='symbol')) %>% filter(!is.na(chromosome))
  # log transform of TPM+1 values 
  ig_expressions = ig_expressions %>% mutate(logTPM=log2(TPM + 1))
  igh_expressions = ig_expressions %>% filter(chain=='IGH')
  igk_expressions = ig_expressions %>% filter(chain=='IGK')
  igl_expressions = ig_expressions %>% filter(chain=='IGL')
  
  pdf(paste0('~/Documents/', sampleID, '_IGH_TPM.pdf'))
  par(mfrow=c(3,1))
  barplot(igh_expressions$logTPM, names=igh_expressions$hgnc_symbol, ylim=c(0, max(ig_expressions$logTPM)), col = rainbow(4)[igh_expressions$VDJC %>% factor(levels=c("V", "D", "J", "C")) %>% as.numeric()], ylab='log2(TPM+1)', main='IGH',  las=2, cex.axis=0.5, cex.names=0.5)
  mtext(sampleID, outer=F, side=3, line=3)
  barplot(igk_expressions$logTPM, names=igk_expressions$hgnc_symbol, ylim=c(0, max(ig_expressions$logTPM)), col = rainbow(4)[igk_expressions$VDJC %>% factor(levels=c("V", "D", "J", "C")) %>% as.numeric()], ylab='log2(TPM+1)', main='IGK',  las=2, cex.axis=0.5, cex.names=0.5)
  barplot(igl_expressions$logTPM, names=igl_expressions$hgnc_symbol, ylim=c(0, max(ig_expressions$logTPM)), col = rainbow(4)[igl_expressions$VDJC %>% factor(levels=c("V", "D", "J", "C")) %>% as.numeric()], ylab='log2(TPM+1)', main='IGL',  las=2, cex.axis=0.5, cex.names=0.5)
  dev.off()
  
}

# get all sample names
sampleNames <- str_replace(basename(Sys.glob("~/Documents/cjyoon/scripts/star/rsem/*.genes.results")), '.genes.results', '')

# plot IG expression graphs
sapply(sampleNames, plot_ig_expressions)


plot_ig_expressions(sampleNames[1])

