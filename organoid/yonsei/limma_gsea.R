# searches anaconda3 fist then other paths for libraries (here due to some fortran error for gplots)
.libPaths( c("/home/users/cjyoon/anaconda3/lib/R/library",  .libPaths()) )
.libPaths()
library(edgeR)
library(limma)
library(Glimma)
library(Homo.sapiens)
library(RColorBrewer)
library(gplots)
library(ggplot2)

# Locate the count files for input
dir = '../../../../Projects/yonsei_organoid/limma/ready_files'
files <- list.files(dir , pattern = "\\.txt")
setwd(dir)

x <- readDGE(files, columns=c(1,3)) 
samplenames <- colnames(x)
group <- as.factor(c('BEM', 'BEM', 'BEM', 'MAT', 'MAT', 'MAT'))
x$samples$group <- group
x$samples

geneid <- rownames(x)
genes <- select(Homo.sapiens, keys=geneid, columns=c("SYMBOL", 'TXCHROM'), keytype="ENTREZID")
# remove those with duplicated annotation
genes <- genes[!duplicated(genes$ENTREZID), ]
# add gene annotation dataframe to the DGEList object
x$genes <- genes
x


# Log transformation into CPM values
# Gene length remain constant across samples, so not necessary to use FPKM. Using CPM 
cpm <- cpm(x)
lcpm <- cpm(x, log=TRUE)

# remove genes those are expressed at low levels
# first, investigate how many are lowly expressed
table(rowSums(x$counts==0)==6)  # genes that are expressed at 0 for all 9 samples


# remove low expressions
keep.exprs <- rowSums(cpm>1) >= 3 
x <- x[keep.exprs, keep.lib.sizes=FALSE]
dim(x)

# Plot raw data
nsamples <- ncol(x)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0, 0.21), las=2, main='', xlab='')
title(main='A. Raw data', xlab='Log-cpm')
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend('topright', samplenames, text.col=col, bty='n')


# now start from x which has low expressed genes filtered out
lcpm <- cpm(x, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0, 0.21), las=2, main='', xlab='')
title(main='B. Filtered data', xlab='Log-cpm')
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
  
}
legend('topright', samplenames, text.col=col, bty='n')


# Normalization
x<- calcNormFactors(x, method='TMM')
x$samples$norm.factors


# MDS plots (PCA)
lcpm <- cpm(x, log=TRUE)
par(mfrow=c(1,2))
col.group <- group
levels(col.group) <- brewer.pal(nlevels(col.group), 'Set1')
col.group <- as.character(col.group)
plotMDS(lcpm, labels=group, col=col.group)
title(main='A. Sample groups')


glMDSPlot(lcpm, labels=group, groups=x$samples[, c(2,4)], launch=F)


# Create Design matrix
design <- model.matrix(~0 + group)
colnames(design) <- gsub('group', '', colnames(design))

# Create Contrast matrix
contr.matrix <- makeContrasts(
  BEMvsMAT = BEM - MAT, 
  levels = colnames(design)
)


# Voom and fit 
v <- voom(x, plot=TRUE)
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit)
summary(decideTests(efit))

# if using a log-fold-changes threshold
tfit <- treat(vfit, lfc=0)
dt <- decideTests(tfit)
summary(dt)


# look at top DE genes
bem_vs_mat <- topTreat(tfit, coef=1, n=Inf)
bem_vs_mat


plotMD(tfit, column=1, status=dt[,1], main=colnames(tfit)[1])
glMDPlot(tfit, coef=1, status=dt, main=colnames(tfit)[1], side.main="ENTREZID", counts=x$counts, group=group, launch=F)


bem_vs_mat.topgenes <- bem_vs_mat$ENTREZID[1:100]
i <- which(v$genes$ENTREZID %in% bem_vs_mat.topgenes)
mycol <- colorpanel(1000, 'blue', 'white', 'red')
pdf('/home/users/cjyoon/Projects/yonsei_organoid/limma/ready_files/top_degenes.pdf', height=14, width=7)
heatmap.2(v$E[i, ], scale='row', 
          labRow=v$genes$SYMBOL[i], labCol=group, 
          col=mycol, trace='none', density.info='none', 
          margin=c(8,6), lhei=c(2,10), dendrogram='column')
dev.off() 
download.file("http://bioinf.wehi.edu.au/software/MSigDB/human_c2_v5p2.rdata", "human_c2_v5p2.rdata", mode = "wb")
download.file('http://bioinf.wehi.edu.au/software/MSigDB/human_H_v5p2.rdata', 'human_h_v5p2.rdata', mode='wb')
download.file('http://bioinf.wehi.edu.au/software/MSigDB/human_c5_v5p2.rdata', 'human_c5_v5p2.rdata', mode='wb')
load("human_c2_v5p2.rdata")
load('human_h_v5p2.rdata')
load('human_c5_v5p2.rdata')
# get indices for Hallmark gene sets
idx_Hallmark <- ids2indices(Hs.H, id=rownames(v))
# get indicies for C2 gene sets
idx_C2 <- ids2indices(Hs.c2, id=rownames(v))
# get indicies for C5 gene sets
idx_C5 <- ids2indices(Hs.c5, id=rownames(v))


cam.bem_vs_mat_Hallmark <- camera(v, idx_Hallmark, design, contrast=contr.matrix[,1])
cam.bem_vs_mat_C2 <- camera(v, idx_C2, design, contrast=contr.matrix[,1])
cam.bem_vs_mat_C5 <- camera(v, idx_C5, design, contrast=contr.matrix[,1])


View(cam.bem_vs_mat_C5)
library(tidyverse)
cam.bem_vs_mat_write_Hallmark = as.tibble(cam.bem_vs_mat_Hallmark) %>% mutate(HALLMARK = row.names(cam.bem_vs_mat_Hallmark)) # added rowname as a separate column for writing
cam.bem_vs_mat_write_C2 = as.tibble(cam.bem_vs_mat_C2) %>% mutate(C2 = row.names(cam.bem_vs_mat_C2)) # added rowname as a separate column for writing
cam.bem_vs_mat_write_C5 = as.tibble(cam.bem_vs_mat_C5) %>% mutate(C5 = row.names(cam.bem_vs_mat_C5)) # added rowname as a separate column for writing

write.table(cam.bem_vs_mat_write_Hallmark, '../../../../Projects/yonsei_organoid/limma/ready_files/gsea/cam_bem_vs_mat_Hallmark_set.txt', sep='\t', quote=F, row.names=F, col.names=T)
write.table(cam.bem_vs_mat_write_C2, '../../../../Projects/yonsei_organoid/limma/ready_files/gsea/cam_bem_vs_mat_C2_set.txt', sep='\t', quote=F, row.names=F, col.names=T)
write.table(cam.bem_vs_mat_write_C5, '../../../../Projects/yonsei_organoid/limma/ready_files/gsea/cam_bem_vs_mat_C5_set.txt', sep='\t', quote=F, row.names=F, col.names=T)

# plot example barcode graph for hallmark
pdf('../../../../Projects/yonsei_organoid/limma/ready_files/gsea/significant_barcodeplot_Hallmark.pdf')
barcodeplot(efit$t[,1], index=idx_Hallmark$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION, index2=idx_Hallmark$HALLMARK_COAGULATION, main="HALLMARK_EPITHELIAL_MESENCHYMAL\n_TRANSITION\nHALLMARK_COAGULATION")
barcodeplot(efit$t[,1], index=idx_Hallmark$HALLMARK_ANDROGEN_RESPONSE, index2=idx_Hallmark$HALLMARK_TGF_BETA_SIGNALING, main="HALLMARK_ANDROGEN_RESPONSE\nHALLMARK_TGF_BETA_SIGNALING")
dev.off()

# plot example barcode graph for C2 gene set
pdf('../../../../Projects/yonsei_organoid/limma/ready_files/gsea/significant_barcodeplot_C2.pdf')
barcodeplot(efit$t[,1], index=idx_C2$GINESTIER_BREAST_CANCER_ZNF217_AMPLIFIED_DN, index2=idx_C2$IGLESIAS_E2F_TARGETS_UP, main="GINESTIER_BREAST_CANCER_ZNF217_AMPLIFIED_DN\nIGLESIAS_E2F_TARGETS_UP")
barcodeplot(efit$t[,1], index=idx_C2$VERHAAK_GLIOBLASTOMA_MESENCHYMAL , index2=idx_C2$YAMAZAKI_TCEB3_TARGETS_UP, main="VERHAAK_GLIOBLASTOMA_MESENCHYMAL\nYAMAZAKI_TCEB3_TARGETS_UP")
dev.off()


# plot barcode graph for C2 gene set in nervous system
pdf('../../../../Projects/yonsei_organoid/limma/ready_files/gsea/significant_barcodeplot_C2_NervousSystem.pdf')
barcodeplot(efit$t[,1], index=idx_C2$REACTOME_SRP_DEPENDENT_COTRANSLATIONAL_PROTEIN_TARGETING_TO_MEMBRANE, index2=idx_C2$REACTOME_NEURONAL_SYSTEM, main="REACTOME_SRP_DEPENDENT_COTRANSLATIONAL_PROTEIN_TARGETING_TO_MEMBRANE\nREACTOME_NEURONAL_SYSTEM")
barcodeplot(efit$t[,1], index=idx_C2$REACTOME_TRANSMISSION_ACROSS_CHEMICAL_SYNAPSES , index2=idx_C2$REACTOME_GABA_SYNTHESIS_RELEASE_REUPTAKE_AND_DEGRADATION, main="REACTOME_TRANSMISSION_ACROSS_CHEMICAL_SYNAPSES\nREACTOME_GABA_SYNTHESIS_RELEASE_REUPTAKE_AND_DEGRADATION")
barcodeplot(efit$t[,1], index=idx_C2$REACTOME_NEUROTRANSMITTER_RELEASE_CYCLE , index2=idx_C2$REACTOME_AMYLOIDS, main="REACTOME_NEUROTRANSMITTER_RELEASE_CYCLE\nREACTOME_AMYLOIDS")
barcodeplot(efit$t[,1], index=idx_C2$REACTOME_NEUROTRANSMITTER_RECEPTOR_BINDING_AND_DOWNSTREAM_TRANSMISSION_IN_THE_POSTSYNAPTIC_CELL , index2=idx_C2$REACTOME_GABA_A_RECEPTOR_ACTIVATION, main="REACTOME_NEUROTRANSMITTER_RECEPTOR_BINDING_AND_DOWNSTREAM_TRANSMISSION_IN_THE_POSTSYNAPTIC_CELL\nREACTOME_GABA_A_RECEPTOR_ACTIVATION")
barcodeplot(efit$t[,1], index=idx_C2$REACTOME_TRAFFICKING_OF_AMPA_RECEPTORS , index2=idx_C2$REACTOME_GLUTAMATE_NEUROTRANSMITTER_RELEASE_CYCLE, main="REACTOME_TRAFFICKING_OF_AMPA_RECEPTORS\nREACTOME_GLUTAMATE_NEUROTRANSMITTER_RELEASE_CYCLE")
dev.off()
