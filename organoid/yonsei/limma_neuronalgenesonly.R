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
dir = '/home/users/cjyoon/Projects/yonsei_organoid/limma/ready_files'
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
print('******summary dt ********')
summary(dt)


# look at top DE genes
bem_vs_mat <- topTreat(tfit, coef=1, n=Inf)
bem_vs_mat

print('****** tfit ******')
print(tfit)

print("****** top table *********")

neuronal_genes = read.table('/home/users/cjyoon/Projects/yonsei_organoid/limma/neuronal_geneonly/neuronal_genes_entrez', sep='\t')
neuronal_genes = neuronal_genes$V1
gene_index = which(tfit$genes$ENTREZID %in% neuronal_genes)
print(gene_index)
de_neuronal_genes = topTable(tfit[gene_index, ], sort.by = 'p', number=Inf, adjust.method="BH")
write.table(de_neuronal_genes, '/home/users/cjyoon/Projects/yonsei_organoid/limma/neuronal_genes_only_de.tsv', sep = '\t', col.names = T, row.names = F,  quote = F)