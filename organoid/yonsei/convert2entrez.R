library(biomaRt)
library(tidyverse)

# define biomart object
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mapping <- as.tibble(getBM(attributes = c("ensembl_gene_id", 'entrezgene'), mart = ensembl))

dir = '/home/users/cjyoon/Projects/yonsei_organoid/limma/'
files <- list.files(dir , pattern = "\\.counts\\.simple\\.txt")
for(file in files){
  count_ensembl = as.tibble(read.table(paste(dir, file, sep=''), sep='\t', header=T))
  joined = left_join(count_ensembl, mapping, by=c('EnsemblGeneID' ='ensembl_gene_id'))
  filtered = joined[!is.na(joined$entrezgene), ]
  output = filtered %>% mutate(EntrezID=entrezgene) %>% dplyr::select(EntrezID, GeneLength, Count)
  output = output[!duplicated(output$EntrezID), ]
  outputpath = paste('/home/users/cjyoon/Projects/yonsei_organoid/limma/', gsub('.counts.simple.txt', '.txt', file), sep='')
  write_delim(output, outputpath, delim='\t', append=F, col_names = T)
}
