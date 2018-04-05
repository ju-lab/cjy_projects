library(tidyverse)
library(stringr)

pcawg_master = read.table('/home/users/cjyoon/Projects/pcawg/pcawg_korea_data_template.txt',  header = T)
pcawg_master = as.tibble(pcawg_master_md5)

# ADD Local MD5sum result to the table
local_aggregate = as.tibble(read.table('/home/users/cjyoon/Projects/pcawg/summarize/localmd5_aggregate.removeconflict.txt', header=F, sep='\t', col.names = c('local_md5sum', 'bam', 'centerName')))
added_local_md5 = (pcawg_master_md5 %>% left_join(local_aggregate, by=c('wgs_alignment_bam_file_name'='bam'))) 

View(added_local_md5 %>% filter(tumor_specimen_count > 1))

# Add TRUTH MD5sum to the table
true_md5 = read.table('/home/users/cjyoon/Projects/pcawg/Pancancer_list_c_TCGA_PSY-1.txt.sourceMD5', header=F, col.names = c('true_md5sum', 'bam'))
true_md5 = as.tibble(true_md5)


added_local_true_md5= added_local_md5 %>% left_join(true_md5, by=c('wgs_alignment_bam_file_name'='bam'))
View(added_local_true_md5)
# Add Sanger Path to the table

sanger_directory = read.table('/home/users/cjyoon/Projects/pcawg/sanger_directory/pcawg_all_abspath_removed_duplicates.txt', header=F, col.names = c('sanger_path'))
sanger_directory = as.tibble(sanger_directory)
sanger_directory_with_basename = sanger_directory %>% mutate(bam = basename(sanger_path))

added_local_true_md5_sangerDirectory = added_local_true_md5 %>% left_join(sanger_directory_with_basename, by=c('wgs_alignment_bam_file_name'='bam'))

# Add column to indicate whether local md5sum matches source (true) md5sum
# If nothing to match then, NA
result = (added_local_true_md5_sangerDirectory %>% mutate(download_status = (local_md5sum == true_md5sum)))
result %>% filter(download_status == FALSE)
result %>% filter(download_status != TRUE)

# Write the result
write_delim(result, '/home/users/cjyoon/Projects/pcawg/summarize/pcawg_korea_download_status_20180326.tsv', delim = '\t')


# status per project
status_per_project = result %>% group_by(dcc_project_code) %>% summarise(downloaded = sum(download_status==TRUE, na.rm=T), failed =sum(download_status==FALSE, na.rm=T), never=sum(is.na(download_status) ))
status_per_project %>% gather(variable, value, downloaded:never) %>% 
ggplot(aes(x=dcc_project_code, y=value, fill = variable)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1))


# complete bams per institute 
result %>% group_by(centerName) %>% summarise(downloaded = sum(download_status==TRUE, na.rm=T))

