# Script to draw CCSER1 deletion in MM samples. 
# Submitted for ASH 2018 abstract
# CCSER1 regions identified from 0~/P/m/a/ccser1$ less ~/reference/human_g1k_v37_starindex/Homo_sapiens.GRCh37.87.gtf | grep CCSER1
# 2018.07.31

library(shape)
library(tidyverse)
CCSER1_exons = tribble(
  ~exon, ~chromosome, ~start, ~end,
  '2', 'chr4', 91229436, 91230759, 
  '3', 'chr4', 91234014, 91234198,
  '4', 'chr4', 91321187, 91321280,
  '5', 'chr4', 91389385, 91389505,
  '6', 'chr4', 91549176, 91549383,
  '7', 'chr4', 91645065, 91645065,
  '8', 'chr4', 91736913, 91736996,
  '9', 'chr4', 91844521, 91844598,
  '10', 'chr4', 92007101, 92007145, 
  '11', 'chr4', 92519723, 92520205
)
CCSER1_exons = CCSER1_exons %>% mutate(median_pos = (start + end)/2)
CCSER1_exons[1,'median_pos'] = CCSER1_exons[2,'median_pos'] -20000
CCSER1_exons[2,'median_pos'] = CCSER1_exons[2,'median_pos'] + 20000
deletions = tribble(
  ~sampleID, ~chromosome, ~start, ~end, 
  '2665', 'chr4', 91978086, 92116128, 
  '3083', 'chr4', 91902510, 92092723, 
  '4593', 'chr4', 91950858, 92009381, 
  '2685', 'chr4', 91485003, 91568663, 
  '4529', 'chr4', 91449459, 91587601, 
)


deletions$ypos = seq(from=1, to=1.9, by=0.2)
par(oma=c(0.1,0.1,0.1,0.1) )
pdf("CCSER1_deletion.pdf")
plot(0, type = "n", xlim=c(min(CCSER1_exons$start)-200000, max(CCSER1_exons$end))+10000, ylim=c(-1, 6), bty='n',xaxt='n',xlab='', yaxt='n', ylab='')
with(CCSER1_exons, rect(start, -0.5, end, 0.5))
with(CCSER1_exons, text(median_pos, -1, exon))
with(deletions, Arrows(start, ypos, end, ypos,  arr.col = "red", code=3, arr.type='triangle'))
text(min(CCSER1_exons$start-100000), 0, 'CCSER1\nExons')
text(min(CCSER1_exons$start-100000), 1.5, 'Deletions')
text(min(CCSER1_exons$start-130000), -1, 'Exon\nnumber')
dev.off()
 #with(deletions, arrows(start, ypos, end, ypos, col='red', angle=30, code=3))
