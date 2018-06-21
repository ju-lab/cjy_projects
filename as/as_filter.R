# This script plots inter variant distance to detect possible gene duplication
# in Ankylosing Spondylitis WGS/WXS projects
# April 6 2018 Jongsoo Yoon (cjyoon@kaist.ac.kr)

library(ggplot2)
library(tidyverse)

# new.txt is the newer version after fixing error. 2018.06.04
data = read.table('/home/users/team_projects/AnkylosingSpondylitis/analysis/freebayes/vcf/annovar_files/B27NE_char.new.txt.fm.txt', sep='\t', header=T, col.names = c('CHROM', 'POS','REF','ALT', 'sample' ,'sample_IDs', 'gene', 'mut_type', 'exac_tot', 'exac_eas', 'exac_max', 'Kova', 'odds_tot', 'odds_eas', 'odds_max', 'odds_kova', 'Father', 'Mother', 'Proband', 'Sibling'))
data = as.tibble(data)

for(chromosome in 1:23){
  print(chromosome)
  mean_dist = data %>% filter(CHROM==chromosome) %>% mutate(BEFOREPOS = lag(POS)) %>% mutate(DIST=abs(POS-BEFOREPOS)) %>% select(CHROM, POS, BEFOREPOS, DIST) %>%
    group_by(CHROM) %>% summarise(mean_dist = mean(DIST, na.rm=T))
  average = as.integer(mean_dist[1, 'mean_dist'])

  filtered_data = data  %>% replace_na(list(odds_max= 1000)) %>% 
    mutate(family_info = Father==TRUE & Mother==FALSE & Sibling==FALSE & Proband==TRUE) %>% 
    filter(CHROM==chromosome, Proband==T) %>% mutate(BEFOREPOS = lag(POS), AFTERPOS=lead(POS)) %>% 
    mutate(BEFOREDIST=abs(POS-BEFOREPOS), AFTERDIST=abs(POS-AFTERPOS)) %>% 
    filter(BEFOREDIST > 100, AFTERDIST>100, family_info==T)
  
  filtered_data %>%
    ggplot(aes(x=POS, y=odds_max, colour=family_info)) + geom_point(aes(size=sample)) + ggtitle(str_c('chromosome: ', as.character(chromosome)))+  scale_y_log10()
  ggsave(str_c('/home/users/team_projects/AnkylosingSpondylitis/Analysis/as_filtering/', chromosome, '_OR.plot.pdf'))
  
  filtered_data_greaterOdds = filtered_data %>% filter(odds_max > 1)
  if(chromosome == 1){
    write_tsv(filtered_data_greaterOdds, append = F, col_names = T, path=str_c('/home/users/team_projects/AnkylosingSpondylitis/Analysis/as_filtering/', chromosome, '_OR_gt1.txt'))
    
  }else{
    write_tsv(filtered_data_greaterOdds, append = F, col_names = F, path=str_c('/home/users/team_projects/AnkylosingSpondylitis/Analysis/as_filtering/', chromosome, '_OR_gt1.txt'))
    
  }
}





