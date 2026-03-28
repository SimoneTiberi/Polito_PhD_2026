rm(list = ls())
#############################################
# load the data:
# 10x droplet-based scRNA-seq PBMC data from 8 Lupus patients 
# before and after 6h-treatment with INF-beta (16 biological samples in total).
#############################################
library(ExperimentHub)
eh = ExperimentHub()
# load data of interest:
sce = eh[["EH2259"]]

sce
colData(sce)
table(colData(sce)$ind, colData(sce)$stim)
# 8 individuals, in both groups
# paired data: before and after stimulation
# 16 biological samples: 8 controls and 8 stimulated

#############################################
# pre-processing:
#############################################
# remove undetected genes
sce = sce[rowSums(counts(sce) > 0) > 0, ]
dim(sce)

# calculate per-cell quality control (QC) metrics
library(scater)
qc = perCellQCMetrics(sce)

# remove cells with few or many detected genes
ol = isOutlier(metric = qc$detected)
mean(ol)
# 6% of the cells removed
sce = sce[, !ol]
dim(sce)

# remove lowly expressed genes
sce = sce[rowSums(counts(sce) > 1) >= 20, ]
dim(sce)

# compute sum-factors & normalize
sce = computeLibraryFactors(sce)
sce = logNormCounts(sce)
sce

#############################################
# Cell type identification:
#############################################
# reference
library(celldex)
reference = celldex::HumanPrimaryCellAtlasData()

# parallel coding on 8 cores
library(BiocParallel)
BPPARAM = MulticoreParam(workers=8)

# cell-type assignment:
library(SingleR)
# it uses log-normalized counts
predicted_cell_type = SingleR(test = sce, 
                     ref = reference, 
                     assay.type.test=1,
                     labels = reference$label.main,
                     BPPARAM = BPPARAM)
# labels indicate the cell-types we aim to identify
table(reference$label.main)

table(predicted_cell_type$labels)
# we add the predictions to the sce object:
colData(sce)$predicted_cell_types = predicted_cell_type$labels

# good agreement between our curated labels and SingleR labels:
table(colData(sce)$predicted_cell_types, colData(sce)$cell)
# Neutrophils and monocytes are both key phagocytes in the innate immune system 
# Megakaryocytes (MKs) and T cells, while distinct in their primary roles (thrombopoiesis vs. adaptive immunity), share surprising similarities, particularly in their ability to act as immune cells

#############################################
# Dimensionality reduction:
#############################################
sce = runUMAP(sce, BPPARAM = BPPARAM)
?runUMAP
# UMAP added to reducedDimNames:
sce

plotUMAP(sce)
# colour by cureted cell type:
plotUMAP(sce, colour_by = "cell")
# colour by individual (8):
plotUMAP(sce, colour_by = "ind")
# turn int into a factor!
colData(sce)$ind = as.factor(colData(sce)$ind)
plotUMAP(sce, colour_by = "ind")
# colour by group:
plotUMAP(sce, colour_by = "stim")

# very clear difference between conditions (and cell types)


# plotTSNE as well (UMAP usually preferred):
# t-SNE: t-stochastic neighbour embedding
# run t-SNE (it take awhile):
# sce = runTSNE(sce, BPPARAM = BPPARAM)
plotTSNE(sce, colour_by = "cell")
plotTSNE(sce, colour_by = "ind")
plotTSNE(sce, colour_by = "stim")

#############################################
# Markers detection:
#############################################
library(scran)
markers = findMarkers(sce, groups = colData(sce)$cell,
                      pval.type="any", BPPARAM = BPPARAM)
# pval.type="any": combined with Simes
markers
# list with 1 element per cell type tested:

head(markers$`B cells`)
# 1 log-FC per cell-type (7 in total)
# 1 overall p-value (aggregated via Simes) and FDR

# marker detection, usually performed:
# - to identify novel markers on known cell types;
# - on unknown clusters (not cell type labeled) to understand their cell type(s)

#############################################
# pseudo-bulk DGE:
#############################################
library(muscat)

# create the biological sample aggregating ind and stim:
sce$id = paste0(sce$stim, sce$ind)
# prepare data for muscat analyses:
sce = prepSCE(sce, 
                kid = "cell", # subpopulation assignments
                gid = "stim",  # group IDs (ctrl/stim)
                sid = "id",   # sample IDs (ctrl/stim.1234)
                drop = TRUE)  # drop all other colData columns

# compute pseudo-bulk counts:
# for DGE, we add (sum) the raw counts (counts):
pb = aggregateData(sce,
                    assay = "counts", fun = "sum",
                    by = c("cluster_id", "sample_id"))
# dim: 5558 16
# 5.5 k genes * 16 samples
pb
# each cell type willl be analyzed separately: we will perform DGE between conditions

pbMDS(pb)

# DGE between conditions for each cell type:
res = pbDS(pb)
# access results table for 1st comparison
tbl = res$table[[1]]
names(tbl)
head(tbl$`B cells`)

# visualize top DE genes:
pbHeatmap(sce, res)

# focus on a single cell-type:
pbHeatmap(sce, res, k = "B cells")

# focus on a single gene (across cell types):
# consistently up-regulated in stim in all cell types:
pbHeatmap(sce, res, g = "ISG20")

# check concordance of top genes across cell type results:

# filter FDR < 5%, abs(logFC) > 2 & sort by adj. p-value
tbl_fil = lapply(tbl, function(u) {
  u = dplyr::filter(u, p_adj.loc < 0.01, abs(logFC) > 2)
  dplyr::arrange(u, p_adj.loc)
})

library(purrr)
de_gs_by_k = map(tbl_fil, "gene")
library(UpSetR)
upset(fromList(de_gs_by_k), 
      nsets = 8)
# CD14+ Monocytes:
# ~half of DGE genes appear only in that cell type
# and ~half in other cell types too
# for the other 7 cell types: most DE genes are in common with other clusters too.

#############################################
# Single-cell DGE:
#############################################
library(distinct)
# Create the design of the study:
  
samples = sce@metadata$experiment_info$sample_id
group = sce@metadata$experiment_info$group_id
design = model.matrix(~group)
# rownames of the design must indicate sample ids:
rownames(design) = samples
design

cpm(sce) = calculateCPM(sce)

# for computational reasons, select 100 genes to analyze:
cpm(sce) = calculateCPM(sce)
set.seed(169612)
sel = sample( rownames(sce), size = 100)
sce_sel = sce[ rownames(sce) %in% sel, ]

# use normalized data, such as counts per million (CPM) or log2-CPM (e.g., logcounts as created via scater::logNormCounts).
set.seed(169612)
library(distinct)
res = distinct_test(x = sce_sel, 
                    name_assays_expression = "logcounts",
                    name_cluster = "cluster_id",
                    name_sample = "sample_id",
                    design = design,
                    column_to_test = 2,
                    min_non_zero_cells = 20,
                    n_cores = 8)
head(top_results(res))

# distinct test is non-parametric
# to facilitate interpretation, we also compute a log2FC of CPMs:
res = log2_FC(res = res,
              x = sce, 
              name_assays_expression = "cpm",
              name_group = "group_id",
              name_cluster = "cluster_id")
head(top_results(res))

# visualize densities:
plot_densities(x = sce,
               gene = "PLSCR1",
               cluster = "B cells",
               name_assays_expression = "logcounts",
               name_cluster = "cluster_id",
               name_sample = "sample_id",
               name_group = "group_id")


plot_densities(x = sce,
               gene = "PLSCR1",
               cluster = "B cells",
               name_assays_expression = "logcounts",
               name_cluster = "cluster_id",
               name_sample = "sample_id",
               name_group = "group_id",
               group_level = TRUE)

plot_cdfs(x = sce,
               gene = "PLSCR1",
               cluster = "B cells",
               name_assays_expression = "logcounts",
               name_cluster = "cluster_id",
               name_sample = "sample_id",
               name_group = "group_id")

# make violin plots:
library(scater)
plotExpression(sce[,sce$ "B cells"],
               features = "PLSCR1", 
               exprs_values = "logcounts",
               log2_values = FALSE,
               x = "sample_id", colour_by = "group_id", ncol = 3) +
  guides(fill = guide_legend(override.aes = list(size = 5, alpha = 1))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#############################################
# Cell type composition differences
#############################################
cell_type_counts = table(sce$cluster_id, sce$sample_id) 
cell_type_counts = unclass(cell_type_counts) 
head(cell_type_counts)
# count data we input: number of cells of each cell type (for every sample)

library(edgeR)
# Attaching some column metadata.
dd = DGEList(cell_type_counts)
dd

# normalize by library size: total number of cells in a sample
dd = calcNormFactors(dd)

# define the model design
design = model.matrix(~ dd$samples$group)
design

# filter rare abundant cell types
keep = filterByExpr(dd, design)
dd = dd[keep,]
table(keep)
# keeping all 8 cell types

# estimate the dispersion:
dd = estimateDisp(dd, design)
# trended dispersion:
head(dd$trended.dispersion)

# fit NB model, keeping fixed the gene dispersion to the estimated values.
fit = glmFit(dd, design)
design
# LRT
# coef specifies the coefficient to test in the design
# in our case: group coefficient
lrt = glmLRT(fit, coef=2)
lrt

res = topTags(lrt, sort.by="PValue", n = Inf)
head(res$table)

# no evidence for any change between in cell type composition

# alternative method, that works with proportions:
library(speckle)
propeller(clusters=sce$cluster_id, sample=sce$sample_id, group=sce$group_id)
# pro: it works and provides actual proportions

# alternative approach: scCODA (not R), based on a hierarchical Dir Multinomial model.
