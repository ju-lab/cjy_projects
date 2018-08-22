library(tidyverse)
library(ggplot2)
library(stringr)

telseq_df = as.tibble(read.table('/home/users/cjyoon/Projects/myeloma/analysis/telseq/result/myeloma_telseq.combined.txt', sep='\t', header=T))

sample_info = as.tibble(read.table('/home/users/cjyoon/Projects/myeloma/sample_info/tumor_match_abs_all.txt', header=F, col.names = c('patientID', 'tumorbam', 'normalbam')))
sample_info_df = sample_info %>% mutate(tumor=str_extract(tumorbam, '[0-9]+')) %>% mutate(normal=str_extract(normalbam, '[0-9]+')) %>% select(c('patientID', 'tumor', 'normal'))

sample_info_tidy = sample_info_df %>% gather('tumor', 'normal', key='sampleType', value='ReadGroup')
sample_info_tidy$patientID = as.character(sample_info_tidy$patientID)
telseq_df$ReadGroup = as.character(telseq_df$ReadGroup)                 

# box plot comparing normal/tumour. Doesn't show which one is connected to which one
telseq_df %>% left_join(sample_info_tidy) %>% select(c('ReadGroup', 'LENGTH_ESTIMATE', 'sampleType', 'patientID')) %>% ggplot(aes(x=sampleType, y=LENGTH_ESTIMATE)) + geom_boxplot() + geom_jitter(width = 0.3) 
ggsave('/home/users/cjyoon/Projects/myeloma/analysis/telseq/result/telseq_t_vs_n.png')

# telomere length comparing tumor and normal, each dot is from same patient.
telseq_df %>% left_join(sample_info_tidy) %>% select(c('LENGTH_ESTIMATE', 'sampleType', 'patientID')) %>% 
  spread(key=sampleType, value=LENGTH_ESTIMATE) %>% 
  ggplot(aes(x=normal, y=tumor)) + geom_point() +
  scale_x_continuous(limits = c(0, 5)) +   scale_y_continuous(limits = c(0, 5))
ggsave('/home/users/cjyoon/Projects/myeloma/analysis/telseq/result/telseq_individuals.png')

telseq_summary= telseq_df %>% left_join(sample_info_tidy) %>% select(c('LENGTH_ESTIMATE', 'sampleType', 'patientID')) %>% 
  spread(key=sampleType, value=LENGTH_ESTIMATE) %>% mutate(ratio = tumor/normal)
t.test(telseq_summary$tumor, telseq_summary$normal, paired = TRUE)

# Cannot say that telomere length has changed from tumor to normal according to Telseq results under paired student's t-test. 