# Normal Breast Cancer Analysis after soft-clip removal 
# and annotation of regions with indel closeness and capture bed regions
# April 4 2018 Jongsoo Yoon (cjyoon@kaist.ac.kr)
library(tidyverse)
library(ggplot2)
data = '/home/users/team_projects/NormalBRCA/analysis/normal_brca_q30_removeclip.count_with_af.indel_bed_annotated.txt'
df = read.table(data, sep='\t', header=T)

df = as.tibble(df)

df = read_delim(data, delim='\t', col_types=cols(close_to_indel=col_logical(), is_in_capture=col_logical()))

# Calculate Mean depth at capture regions
df %>% filter(is_in_capture==T) %>% group_by(bam) %>% summarise(mean_depth=mean(total))
# # A tibble: 2 x 2
# bam                                                      mean_depth
# <chr>                                                         <dbl>
#   1 114846-D-19.sorted.whole.postdedup.rg.q30.removeclip.bam      1938.
# 2 92247-D-29.sorted.whole.postdedup.rg.q30.removeclip.bam       3207.

base_pairing<-function(base){
  switch(base, 
         A='T', 
         T='A',
         G='C', 
         C='G')
}
base_change_tibble= tribble(
  ~basechange, ~base_change_std, 
  "A>T", "T>A", 
  "A>G", "T>C", 
  "A>C", "T>G",
  "A>A", "T>T",
  "A>N", "T>N",
  "T>A", "T>A",
  "T>G", "T>G", 
  "T>C", "T>C", 
  "T>T", "T>T",
  "T>N", "T>N",
  "G>C", "C>G", 
  "G>T", "C>A", 
  "G>A", "C>T",
  "G>G", "G>G", 
  "G>N", "C>N",
  "C>G", "C>G", 
  "C>T", "C>T", 
  "C>A", "C>A", 
  "C>C", "C>C", 
  "C>N", "C>N"
)


# base_change_standardize <- function(row){
#   if (str_detect(row[, 'basechangeString'], '^[G|C]')){
#     return(row$basechangeString)
#   }else{
#     splitresult = unlist(str_split(row[, 'basechangeString, '>')) 
#     ref = base_pairing(splitresult[1])
#     alt = base_pairing(splitresult[2])
#     return(str_c(ref, '>',alt))
#   }
}

df_withbc =(df %>% mutate(basechange = str_c(reference_allele, '>', minor_allele)))
df_withbc_std = (df_withbc %>% left_join(base_change_tibble))
df_withbc_std %>% filter(close_to_indel==F & is_in_capture==T & minor_allele!='N' & chromosome=='chr1' & total > 1000) %>% 
  ggplot(aes(x=position, y=minor_af, color=base_change_std)) + geom_point()

# Calculate background error rate and std
background_err_rate = df_withbc_std %>% filter(minor_af < 0.4) %>% group_by(bam) %>% summarize(mean = mean(minor_af), std = sd(minor_af))

df_withbc_std_zscore = df_withbc_std %>% left_join(background_err_rate) %>% mutate(zscore = (minor_af-mean)/std) 

df_withbc_std_zscore %>% filter(close_to_indel==F & is_in_capture==T & minor_allele!='N' & chromosome=='chr1' & minor_af < 0.4) %>% 
  ggplot(aes(x=position, y=zscore, color=base_change_std)) + geom_point()
