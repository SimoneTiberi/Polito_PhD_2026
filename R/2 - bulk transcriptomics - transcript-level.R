#############################################
# DAS - Dirichlet Multinomial approach
#############################################
rm(list = ls())
library(BANDITS)
?BANDITS

# load gene_to_transcript matching:
data("gene_tr_id", package = "BANDITS")
head(gene_tr_id)

# We define the design of the study
samples_design = data.frame(sample_id = paste0("sample", seq_len(4)),
                            group = c("A", "A", "B", "B"))
samples_design
# 4 samples that belong to 2 groups

# load the pre-computed data:
data("input_data", package = "BANDITS")
input_data
# small example dataset with 40 genes

# Filter lowly abundant genes:
input_data = filter_genes(input_data, min_counts_per_gene = 20)

# input, not a count matrix, but equivalence classes and their counts
# this gene has 3 ECs (columns) and 2 transcripts (rows):
input_data@classes[[3]]
# columns = EC, rows = transcripts (1 if a transcript is present in the EC)
# 1st EC has both trasncripts, 2nd EC only ENST00000378891, 3rd EC only ENST00000378888

# EC counts: 
input_data@counts[[3]]
# columns = samples; rows = ECs

# load the pre-computed precision estimates:
data(precision, package = "BANDITS")
plot_precision(precision)
# precision/dispersion parameter estimates
# used to formulate an informative prior (normal) for log-precision

# fit the Dirichlet-Multinomial model
# and Test for DAS via a Wald-test
# for each gene, we run an MCMC
set.seed(169612)
results = test_DTU(BANDITS_data = input_data,
                   precision = precision$prior,
                   samples_design = samples_design,
                   R = 10^4, burn_in = 2*10^3, 
                   n_cores = 2,
                   gene_to_transcript = gene_tr_id)
results

# most significant genes:
head(top_genes(results))
# most significant transcripts:
head(top_transcripts(results))

# 2-level test: genes and individual transcripts

# study the most significant gene:
top_gene = top_genes(results, n = 1)
gene(results, top_gene$Gene_id)
# first 2 transcripts are changing between conditions (see transcript_results)

# plot the transcripts relative abundance (pi) in each group
plot_proportions(results, top_gene$Gene_id)
# this gene expresses 4 transcripts
# group B usually expresses the 1st transcript (75% of the times)
# while group A usually expresses the 2nd transcript (70% of the times)

# CI = Wald-type CI based on the sd estimated from the posterior chains.

#############################################
# Load and filter data (same data used in 1st lab):
#############################################
rm(list = ls())
load("data/transcript_counts.RData")
ls()
dim(Tr_counts)
# 207 k transcripts and 17 samples

# we remove the 2 outliers we identified: A_12 and B_11
sel_out = which(samples %in% c("A_12", "B_11"))

samples = samples[-sel_out]
group = group[-sel_out]
Tr_counts = Tr_counts[, -sel_out]

# filter lowly abundant transcripts based on samples kept:
sel_tr = rowSums(Tr_counts) >= 20
mean(sel_tr)
Tr_counts = Tr_counts[ sel_tr, ]
dim(Tr_counts)
# 68 k transcripts

head(Tr_counts)
# ENST00000415118 ... transcript ids
# what else do we need?

# gene ids (transcript-gene match) - what transcript belong to each gene
# stored here:
load("data/tx2gene.RData")
head(tx2gene)

matches = match( rownames(Tr_counts), tx2gene$transcript)
# double check match is correct:
head(rownames(Tr_counts)); head(tx2gene$transcript[matches])
gene_id = as.character( tx2gene$gene[matches] )
# NA -> not in our tx2gene database

# remove NAs in matches:
Tr_counts = Tr_counts[!is.na(matches),]
gene_id = gene_id[!is.na(matches)]

# we can remove single-isoform genes:
# we cannot study difference in alternative splicing, because it does not appear in single-isoform genes!
tab = table(gene_id)
head(tab)
multi_iso_genes = names(tab)[tab > 1.5]

# filter multi-isoform genes in our data:
sel_genes = which(gene_id %in% multi_iso_genes)
Tr_counts = Tr_counts[ sel_genes ,]
gene_id = gene_id[ sel_genes ]
dim(Tr_counts); length(gene_id)
# 57 k isoforms
length(unique(gene_id))
# 12 k genes
rm(sel_genes)

# keep Genes with >= 20 counts per condition (across all samples):
# we aim to estimate proportions in each condition
Tr_counts = as.data.frame(Tr_counts)
counts_split = split(Tr_counts, gene_id)
# transcript counts for each gene
head(counts_split, 2)
sel = sapply(counts_split, function(x){
  # sel samples in group A:
  sel_a = group == "A"
  # TRUE only if min 10 counts in each group (A and B)
  (sum(x[sel_a]) >= 10) & (sum(x[-sel_a]) >= 10)
})
mean(sel)
# 0.988361
genes_kept = unique(names(sel[sel]))
length(genes_kept)
# 12 k genes
sel_genes = gene_id %in% genes_kept
Tr_counts = Tr_counts[ sel_genes ,]
gene_id = gene_id[ sel_genes ]

dim(Tr_counts); length(gene_id)
# 57 k isoforms
length(unique(gene_id))
# 12 k genes

#############################################
# DAS - DEXSeq
#############################################
library(DEXSeq)
set.seed(169612)

# Load truth table:
design = data.frame( row.names = samples,
                     condition = group )
design
# we can add further covariates in the design if needed 
# (e.g., confounders, batch effects, important variables, etc...)

# create DEXSeq data object
dxd = DEXSeqDataSet(countData = round( Tr_counts ),
                    sampleData = design,
                    design= ~ sample + exon + condition:exon,
                    featureID = rownames(Tr_counts),
                    groupID = gene_id, 
                    transcripts = rownames(Tr_counts))

# parallel coding - use 4 cores:
library(BiocParallel)
BPPARAM = MulticoreParam(workers=8)

# estimate normalization constants:
dxd = estimateSizeFactors(dxd)
# shrinkage dispersion estimation:
dxd = estimateDispersions(dxd, BPPARAM = BPPARAM)
# test for DIU (or DEU):
dxd = testForDEU(dxd, reducedModel = ~sample + exon, BPPARAM = BPPARAM)
# extract results
res = DEXSeqResults(dxd, independentFiltering = FALSE)
# limitations:
# i) 1 result per transcript
# ii) no transcript relative abundance
head(res)
# aggregated to gene-level (min p-val + correction):
qval = perGeneQValue(res)
head(qval)
res_DEXSeq = data.frame(gene = names(qval), FDR = qval)

hist(res_DEXSeq$FDR)
# does it look odd?

# we know the distribution of p-vals under H0,
# so looking at p-vals, we expect a flat hist on the right (towards 1)
# but adjusted p-vals have no theorethical distribution (under H0 or H1)
#############################################
# Over-Representation Analysis (ORA) - DAS
#############################################
# choose a threshold (FDR < 0.01) and select significant genes:
genes = res_DEXSeq$gene[res_DEXSeq$FDR < 0.01]
length(genes)
# 2.4 k significant genes
head(genes)

# ORA analysis (exact hyper-geometric over-representation test)
library(clusterProfiler)
EGO_DEXSeq = enrichGO(gene         = genes,
               OrgDb        = 'org.Hs.eg.db',
               keyType      = 'ENSEMBL',
               pvalueCutoff = 1,
               ont = "ALL",
               pool = TRUE)
# ont specifies what onthologies to consider for the gene sets

# order significant pathways first:
EGO_DEXSeq@result = EGO_DEXSeq@result[order(EGO_DEXSeq@result$p.adjust), ]
dim(EGO_DEXSeq@result)
# 8 k pathways tested

# significant pathways:
sum(EGO_DEXSeq@result$p.adjust < 0.01)
# ~1 k significant pathways (FDR < 0.01)

head(EGO_DEXSeq@result[,-which(colnames(EGO_DEXSeq@result) == "geneID")])

library(DOSE)
dotplot(EGO_DEXSeq, showCategory=10, x = "p.adjust")
# functional pathways (several binding terms)

save(res_DEXSeq, EGO_DEXSeq, file = "output/DEXSeq.RData")

# GSEA usually rank by -log(p) * sign(FC)
# but FC not defined in DAS.

#############################################
# DGE and DAS overlap 
#############################################
rm(list = ls())
load("output/edgeR.RData")
load("output/DEXSeq.RData")
ls()

#keep gene and FDR only for edgeR results:
res_edgeR = data.frame(gene = rownames(res$table), FDR_edgeR = res$table$FDR)
colnames(res_DEXSeq)[2] = "FDR_DEXSeq"
# merge edgeR and DEXSeq results:
# all = FALSE - keep intersection only
DF = merge(res_edgeR, res_DEXSeq, by = "gene", all = FALSE)
head(DF)

# check if edgeR and DEXSeq results are associated:
plot(DF$FDR_edgeR, DF$FDR_DEXSeq)
cor(DF$FDR_edgeR, DF$FDR_DEXSeq)
# -0.08503758

# check if significant results are similar:
# we need to dicotomize FDR (> vs. < 0.01):
# UpSetR = Venn diagram but easier to interpret with many methods
DF_upset = data.frame(edgeR = ifelse(DF$FDR_edgeR < 0.01, 1,0),
                      DEXSeq = ifelse(DF$FDR_DEXSeq < 0.01, 1,0) )
library(UpSetR)
upset( DF_upset )
# many more DGE than DAS
# among DAS: 46% also are DGE

# a combination of 2 reasons:
# 1) some "key" genes display both changes;
# 2) more counts = more statistical power to detect both.

# in this case: a lot of overlap, but in general
# results often orthogonal: DGE and DAS are distinct biological processes

#############################################
# ORA-DGE and ORA-DAS overlap
#############################################
rm(list = ls())
load("output/edgeR.RData")
load("output/DEXSeq.RData")
ls()

#keep gene and FDR only for edgeR results:
res_edgeR  = data.frame(GO_ID = EGO@result$ID, 
                        FDR_edgeR = EGO@result$p.adjust)
res_DEXSeq = data.frame(GO_ID = EGO_DEXSeq@result$ID, 
                        FDR_DEXSeq = EGO_DEXSeq@result$p.adjust)
# merge edgeR and DEXSeq results:
# all = FALSE - keep intersection only
DF = merge(res_edgeR, res_DEXSeq, by = "GO_ID", all = FALSE)
head(DF)

# check if edgeR and DEXSeq results are associated:
plot(DF$FDR_edgeR, DF$FDR_DEXSeq)
cor(DF$FDR_edgeR, DF$FDR_DEXSeq)
# 0.5077607

# check if significant results are similar:
DF_upset = data.frame(edgeR = ifelse(DF$FDR_edgeR < 0.01, 1,0),
                      DEXSeq = ifelse(DF$FDR_DEXSeq < 0.01, 1,0) )
library(UpSetR)
upset( DF_upset )
# a lot more overlap at pathway level:
# ~ 1/3 of pathways display DGE only
# ~ 1/3 of pathways display DAS only
# ~ 1/3 of pathways display DGE + DAS

# ~half of DAS genes, also show DGE
# ~half of DGE genes, also show DAS

# considering groups of genes with similar functions, we have more similar results
