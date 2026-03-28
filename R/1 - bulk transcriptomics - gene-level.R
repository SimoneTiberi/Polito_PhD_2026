#############################################
# load the data:
#############################################
# here, we use a filtered dataset.
rm(list = ls())
load("data/gene_counts.RData")
ls()
dim(Gene_counts)
# 36 k genes (including NON-protein coding genes)
# 17 samples
# humans ~ 20k protein-coding genes

head(Gene_counts)
# actual counts are discrete, but these are continuous.
# because they are estimates, due to multi-mapping reads
# these are gene counts: multi-mapping across genes (less common than across transcripts)

samples
group
table(group)
# we aim to compare these 2 groups

# define a colour for each group"
cols = ifelse(group == "A", "blue", "red")
cols

# Inspect gene counts:
hist(Gene_counts)
mean(Gene_counts == 0)
# 0.4486605

# filter lowly abundant genes:
# keep genes with min 20 counts across all samples:
sel_gene = rowSums(Gene_counts) >= 20
mean(sel_gene)
Gene_counts = Gene_counts[ sel_gene, ]
Gene_length = Gene_length[ sel_gene, ]

dim(Gene_counts)
# ~22 k  genes

#############################################
# Exploratory plots
#############################################
library(ggplot2)

# normalize Transcript_counts:
library(edgeR);
dd = DGEList(counts = Gene_counts, samples = samples)
dd = calcNormFactors(dd)
dd$samples
# lib.size (total counts in sample; M_i in slides); same as:
colSums(Gene_counts)
# difference close to 0:
colSums(Gene_counts) - dd$samples$lib.size

# MDS plot (similar to PCA, but non-linear dimensionality reduction)
m3 = plotMDS(dd, labels = group, 
             col = cols, main = "MDS")
# 1) clear separation between groups
# 2) 2 outliers
m3 = plotMDS(dd, labels = samples, 
             col = cols, main = "MDS")
# outliers are samples B_11 and A_12

# PCA plot:
# we do not normalize the data (0 mean and 1 sd)
# because we do NOT want all covariates (genes) to have the same influence!
# usually, we work on log-CPMs: the log decreases the importance of the highly abundant genes
library(ggfortify)
# normalize Transcript_counts (log-CPMs)
logcounts = cpm(dd, log = TRUE)
colnames(logcounts) = samples
pca = prcomp(t(logcounts))
rownames(pca$x) = samples
autoplot( pca, col = cols, 
          label = TRUE) + theme_bw()
# different dimentionality reduction but similar conclusion.
# PCA VERY effective:
# ~20 k covariates (genes) -> top 2 components capture ~2/3 of the variance
# PC1 -> mainly separates outliers
# PC2 -> mainly separates groups

# k-means clustering:
set.seed(169612)
km = kmeans(t(logcounts), 2, 
            nstart = 100, iter.max = 1000)
# compare kmeans groups with clusters:
table(km$cluster, group)
# group 1: outliers (samples B_11 and A_12); and group 2: all remaining samples
km$cluster 

set.seed(169612)
km = kmeans(t(logcounts), 3,
            nstart = 100, iter.max = 1000)
table(km$cluster, group)
# 2 groups for outliers
# 1 group for A (excluding A_12)
# 1 group for B (excluding B_11)

# same conclusions: 2 groups are clearly distinct,
# and there are 2 outliers

# heatmap based on the most variable genes (across samples):
# HVGs -> most informative in discerning samples
var_genes = apply(logcounts, 1, var)
# we select the 10^3 most variable log-cpm
select_var = names(sort(var_genes, decreasing=TRUE))[1:1000]
# Subset logcounts matrix (only keep the most variables log-cpm)
highly_variable_lcpm = logcounts[select_var,]
# we apply a hierarchical clustering on the samples and on the rows.
colnames(highly_variable_lcpm) = samples

library(corrplot); library( gplots) 
heatmap.2(highly_variable_lcpm, trace="none", 
          main="Top 1000 variable Genes")
# again same conclusions as above.
# 2 outliers have a very distinct expression profile.

# Correlation plots between samples's expression patterns:
logcounts = logcounts[,order(group)]
corrplot( cor(logcounts), method="color" )
# ordering by group helps to see group-level structures
# higher correlation INSIDE A and B groups
# also A samples are more homogeneous (less sample to sample variability)
# compared to B samples

#############################################
# Filter data:
#############################################
# we remove the 2 outliers: A_12 and B_11
sel_out = which(samples %in% c("A_12", "B_11"))

samples = samples[-sel_out]
group = group[-sel_out]
Gene_counts = Gene_counts[, -sel_out]
Gene_length = Gene_length[, -sel_out]

# re-filter lowly abundant genes based on samples kept:
sel_gene = rowSums(Gene_counts) >= 20
mean(sel_gene)
Gene_counts = Gene_counts[ sel_gene, ]
Gene_length = Gene_length[ sel_gene, ]
dim(Gene_counts)
# 21 k genes

#############################################
# DGE - compute scaling factors - explained in:
# Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences
# Charlotte Soneson, Michael I. Love, and Mark D. Robinson
#############################################
# gene length is a matrix: gene length should not change with sample id: why?
head(Gene_length)
# because it is an average of transcript lengths, weighted by the transcript proportions in each sample.

# Obtaining per-observation scaling factors for length, adjusted to avoid
# changing the magnitude of the counts.
normMat = Gene_length/exp(rowMeans(log(Gene_length)))
normCts = Gene_counts/normMat

# Computing effective library sizes from scaled counts, to account for
# composition biases between samples.
library(edgeR)
set.seed(169612)
eff.lib = calcNormFactors(normCts) * colSums(normCts)

# Combining effective library sizes with the length factors, and calculating
# offsets for a log-link GLM.
normMat = sweep(normMat, 2, eff.lib, "*")
normMat = log(normMat)

#############################################
# DGE - analysis
#############################################
# use original NOT-NORMALIZE COUNTS
y = DGEList(Gene_counts, group = group)
# include normalization factor
y = scaleOffset(y, normMat)
# design with library size 
y$samples
# actual counts:
y$counts
# normalization factor:
y$offset

# design
design = model.matrix(~ group)
design # design matrix of our model, the only covariate is the group (A or B)
# filter lowly abundant genes (if not enough Gene_counts: results are not reliable)
?filterByExpr
# several min to be respected:
keep = filterByExpr(y, design)
mean(keep); sum(keep)
# we keep ~12.5 k genes
y = y[keep, ]
dim(y)

# estimate the dispersion:
y = estimateDisp(y, design)
# trended dispersion:
head(y$trended.dispersion)

# fit NB model, keeping fixed the gene dispersion to the estimated values.
fit = glmFit(y, design)
design
# LRT
# coef specifies the coefficient to test in the design
# in our case: group coefficient
lrt = glmLRT(fit, coef=2)
lrt

# extract results, sorted by significance
res = topTags(lrt, sort.by="PValue", n = Inf)
head(res$table)
# logFC = log2( FC ) between groups
# log simplifies FC interpretation:
# logCPM (how abundant the gene it)
# LRT statistic
# LRT p-value
# FDR = BH adjusted p-values

# plot raw p-values
hist(res$table$PValue)
# clear separation between groups: lots of significant results

table(res$table$FDR < 0.01)
# ~5.8 k significant genes at the 1% FDR level!
# of these results, I expect ~58 (1%) to be FPs

# summarize log2-FC and p-value in a volcano plot
# significance and magnitude:
library(EnhancedVolcano)
EnhancedVolcano(res$table,
                lab = rownames(res$table),
                x = 'logFC',
                y = 'PValue')
?EnhancedVolcano

# we can change p-value and logFC cutoff:
EnhancedVolcano(res$table,
                lab = rownames(res$table),
                x = 'logFC',
                y = 'PValue',
                pCutoff = 10e-6, #def: 10e-5
                FCcutoff = 2) # def: 1

#############################################
# Over-Representation Analysis (ORA) - DGE
#############################################
# choose a threshold (FDR < 0.01) and select significant genes:
genes = rownames(res$table)[res$table$FDR<0.01]
length(genes)
head(genes)

# ORA analysis (exact hyper-geometric over-representation test)
library(clusterProfiler)
EGO = enrichGO(gene         = genes,
               OrgDb        = 'org.Hs.eg.db',
               keyType      = 'ENSEMBL',
               pvalueCutoff = 1,
               ont = "ALL",
               pool = TRUE)
# ont specifies what onthologies to consider for the gene sets

# order significant pathways first:
EGO@result = EGO@result[order(EGO@result$p.adjust), ]
dim(EGO@result)
# 8.7 k pathways tested (instead of ~12.5 k genes)

# significant pathways:
sum(EGO@result$p.adjust < 0.01)
# 1003 significant pathways (FDR < 0.01)
# 5.8 k  significant genes

head(EGO@result[,-which(colnames(EGO@result) == "geneID")])
# Description: knowledge of the pathway
# GeneRatio: significant genes in pathway/all significant genes
# k/K (127/5570)
# BgRatio: all genes in pathway/all genes in DB (23472)
# n/N (183/23472)
# Pathway 1: has 183 genes (n), of which 127 significant (k).
# in total we have 23472 genes (N), of which 5570 significant (K).
# FoldEnrichment: (k/n) / (K/N)
# fraction of significant genes in pathways (0.69 in Pathway 1)
# over fraction of significant genes in database (0.24)
# Pathway 1: (127/183) / (5570/23472) = 0.69 / 0.24 = 2.9
# Pathway 1 has 2.9 times more significant genes than expected under H0.
# p-value: Hypergeometric test p-value
# p.adjust: BH FDR correction
# qvalue: Storey’s FDR correction

library(DOSE)
dotplot(EGO, showCategory=10, x = "p.adjust")
# functional pathways (several binding terms)

save(res, EGO, file = "output/edgeR.RData")

#############################################
# fgsea GSEA analysis
#############################################
library(fgsea)
# fgseaRes will use the score to order results
score = -log10(res$table$PValue) * sign(res$table$logFC)
hist(score)
names(score) = rownames(res$table)
# centre = non-DGE genes
# right: up-regulated genes
# left: down-regulated genes
# as shown by the volcano plot: more down than up

# download reference pathways:
library(msigdbr)
pathwaysDF = msigdbr(species = "human")
pathwaysDF[1,]
# key: ensembl_gene and  gs_name:
# for each pathway (gs_name) list all its genes:
pathways = split(pathwaysDF$ensembl_gene, pathwaysDF$gs_name)
pathways[1]

length(pathways)
# 35 k pathways: issue -> largely overlapping

# GSEA (it takes a few mins):
fgseaRes = fgsea(pathways = pathways, 
                 stats    = score,
                 eps = 0)
# sort results by significance:
fgseaRes = fgseaRes[order(pval), ]

# hist of p-vals:
hist(fgseaRes$pval)

# visualize top results:
head(fgseaRes)

# NES = normalized enrichment score
# for positive ES: NES = ES / average_positive_ES
# for negative ES: NES = ES / average_negative_ES

# store top 10 UP and DOWN regulated pathways
top_UP = fgseaRes$pathway[fgseaRes$ES > 0][1:10]
top_DOWN = fgseaRes$pathway[fgseaRes$ES < 0][1:10]

# plot enrichment of specific pathways
# plot all genes analyzed
# lines = genes in pathway
# mostly to left - pathway UP:
plotEnrichment(pathways[[top_UP[1]]],
               score)

# mostly to right - pathway DOWN:
plotEnrichment(pathways[[top_DOWN[1]]],
               score)

# plot multiple pathways:
plotGseaTable(pathways[top_UP],
              score, fgseaRes)

plotGseaTable(pathways[top_DOWN],
              score, fgseaRes)
