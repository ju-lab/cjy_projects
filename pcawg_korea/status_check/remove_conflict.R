library(tidyverse)
data = read_delim('localmd5_aggregate.txt', delim='\t', col_names=c('local_md5sum', 'bam', 'centerName'))

two_distinct_md5 = data %>% distinct(local_md5sum, bam) %>% group_by(bam) %>% summarise(count = n()) %>% filter(count != 1)
distinct_results = data %>% distinct(local_md5sum, bam, centerName) %>% filter(!(bam %in% two_distinct_md5$bam))
write_delim(distinct_results, 'localmd5_aggregate.removeconflict.txt', delim='\t', col_names=T)

