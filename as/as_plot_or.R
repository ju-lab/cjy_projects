# Plots per chromosome ODDS Ratio of observing a variant for AS. 
###########################################################################
# README
# WILL ADD INFO ON FAMILY TRANSMISSION
# 
###########################################################################
# March 30 2018
# Jongsoo Yoon


library(tidyverse)
data = read.table('/home/users/team_projects/AnkylosingSpondylitis/analysis/freebayes/vcf/annovar_files/B27NE_char.txt.fm.txt', header=T, col.names = c('CHROM', 'POS','REF','ALT', 'sample' ,'sample_IDs', 'gene', 'exac_tot', 'exac_eas', 'exac_max', 'Kova', 'odds_tot', 'odds_eas', 'odds_max', 'odds_kova', 'Father', 'Mother', 'Proband', 'Sibling'))
data = as.tibble(data)
data
data %>% arrange(desc(odds_max))

data %>% mutate(family_info = Father==TRUE & Mother==FALSE & Sibling==FALSE & Proband==TRUE)
for(chrom in 1:23){
  print(chrom)
  data %>% mutate(family_info = Father==TRUE & Mother==FALSE & Sibling==FALSE & Proband==TRUE) %>% filter(CHROM==chrom) %>% ggplot(aes(x=POS, y=odds_max, colour=family_info)) + geom_point(aes(size=sample), alpha=0.3) + ggtitle(str_c('chromosome: ', as.character(chrom))) 
  print(str_c('/home/users/team_projects/AnkylosingSpondylitis/odd_ratio_chrom', as.character(chrom), '.png'))
                    ggsave(str_c('/home/users/team_projects/AnkylosingSpondylitis/odd_ratio_chrom', as.character(chrom), '.png'))
}

chrom = 12
data %>% filter(CHROM==chrom, Mother==F, Sibling==F, Father==T) %>% ggplot(aes(x=POS, y=odds_max, colour=Father)) + geom_point(aes(size=sample)) + ggtitle(str_c('chromosome: ', as.character(chrom))) 

data %>% mutate(family_info = Father==TRUE & Mother==FALSE & Sibling==FALSE & Proband==TRUE) %>% filter(CHROM==chrom) %>% ggplot(aes(x=POS, y=odds_max, colour=family_info)) + geom_point(aes(size=sample), alpha=0.3) + ggtitle(str_c('chromosome: ', as.character(chrom))) + theme_bw()

