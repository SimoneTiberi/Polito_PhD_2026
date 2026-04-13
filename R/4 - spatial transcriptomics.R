rm(list = ls())

library(DESpace)
library(ggplot2)
library(ggforce)
library(SpatialExperiment)
library(spatialLIBD)
library(ExperimentHub)
library(scuttle)

################ ################ ################ ################ 
# load data and packages
################ ################ ################ ################ 
# Connect to ExperimentHub
ehub = ExperimentHub::ExperimentHub()
# Download the full real data (about 2.1 GB in RAM) use:
spe_all = spatialLIBD::fetch_data(type = "spe", eh = ehub)
# Specify column names of spatial coordinates in colData(spe) 
coordinates = c("array_row", "array_col")
# Specify column names of spatial clusters in colData(spe) 
# manually curated clusters:
cluster_col = 'layer_guess_reordered'

# the dataset has multiple samples, extract one only:
spe = spe_all[, colData(spe_all)$sample_id == '151673']
rm(spe_all)

# Remove spots missing annotations
spe = spe[, !is.na(spe[[cluster_col]])]

################ ################ ################ ################ 
# remove low quality spots and low abundant genes:
################ ################ ################ ################ 
spe = scuttle::addPerCellQC(spe,)

# Remove combined set of low-quality spots:
# min 25 counts across all genes per spot
# min 25 non-zero genes detected per spot
# expr_chrM_ratio < 0.3
# number of cells estimated in each spot <= 15
spe = spe[, !(colData(spe)$sum < 25 | # spot-library size
                colData(spe)$detected < 25 | # number of expressed genes
                colData(spe)$expr_chrM_ratio > 0.3| # mitochondrial expression ratio
                colData(spe)$cell_count > 15)] # number of cells per spot
# Then, we discard lowly abundant genes, which were detected in less than 20 spots.
# from 3611 to 3594 spots


# Remove lowly abundant genes, with less than 20 non-zero spots:
qc_low_gene = rowSums(assays(spe)$counts > 0) >= 20
spe = spe[qc_low_gene,]
spe
# from 33538 to 15110 genes

# visualize manually annotated spatial domains:
CD = as.data.frame(colData(spe))
ggplot(CD, 
       aes(x=array_col,y=array_row, 
           color=factor(layer_guess_reordered))) +
  geom_point() + 
  theme_void() + scale_y_reverse() + 
  theme(legend.position="bottom") + 
  labs(color = "", title = paste0("Manually annotated spatial clusters"))

################ ################ ################ ################ 
# Compute Spatial clustering via BayesSpace
################ ################ ################ ################ 
if(FALSE){
  library(BayesSpace)
  library(dplyr)
  
  set.seed(101)
  dec = scran::modelGeneVar(spe)
  top = scran::getTopHVGs(dec, n = 2000)
  
  set.seed(102)
  spe = scater::runPCA(spe, subset_row=top)
  
  ## Add BayesSpace metadata
  spe = spatialPreprocess(spe, platform="Visium", skip.PCA=TRUE)
  
  q = 7  # Number of clusters
  d = 15  # Number of PCs
  
  ## Run BayesSpace clustering
  set.seed(104)
  spe = spatialCluster(spe, q=q, d=d, platform='Visium',
                       nrep=50000, gamma=3, save.chain=TRUE)
  
  # save(spe, file = "output/spe.RData")
  # load("output/spe.RData")
  
  CD = as.data.frame(colData(spe))
  ggplot(CD, 
         aes(x=array_col,y=array_row, 
             color=factor(spatial.cluster))) +
    geom_point() + 
    theme_void() + scale_y_reverse() + 
    theme(legend.position="bottom") + 
    labs(color = "", title = paste0("BayesSpace clusters"))
}

################ ################ ################ ################ 
# Identify Spatially Variable Genes via DESpace
################ ################ ################ ################ 
library(DESpace)
set.seed(169612)
results = svg_test(spe = spe,
                   cluster_col = cluster_col, 
                   verbose = TRUE)
head(results$gene_results, 3)

# plot 6 top SVGs:
feature = results$gene_results$gene_id[seq_len(6)]

FeaturePlot(spe, feature, 
            coordinates = coordinates, 
            ncol = 3, title = TRUE)
# add cluster boundaries:
FeaturePlot(spe, feature, 
            coordinates = coordinates, 
            annotation_cluster = TRUE,
            cluster_col = cluster_col, 
            cluster = 'all', title = TRUE)

# Individual region test:
set.seed(169612)
cluster_results = individual_svg(spe, 
                                 edgeR_y = results$estimated_y,
                                 cluster_col = cluster_col)
names(cluster_results)
class(cluster_results)
# list: 1 element per spatial cluster
# logFC (log2-FC) compares expression in WM vs. non-WM:
head(cluster_results$WM)

# merge global tissue and individual cluster results (FDR and logFC):
merge_res = top_results(results$gene_results, cluster_results)
head(merge_res,3)

# add individual cluster FDR only (not logFC):
merge_res = top_results(results$gene_results, cluster_results, 
                        select = "FDR")
head(merge_res,3)

# Check top genes for a specific cluster: WM
results_WM_both = top_results(cluster_results = cluster_results, 
                              cluster = "WM", 
                              high_low = "both")
# high_low = "both" -> creates 2 lists, 1 for high abundance in WM, and 1 for low abundance in WM:
str(results_WM_both)
head(results_WM_both$high_genes, 3)
head(results_WM_both$low_genes, 3)

# Plot SVGs with higher than average abundance in WM
feature = rownames(results_WM_both$high_genes)[seq_len(3)]
FeaturePlot(spe, feature, cluster_col = cluster_col, 
            coordinates = coordinates, cluster = 'WM', 
            legend_cluster = TRUE, annotation_cluster = TRUE, 
            linewidth = 0.6, title = TRUE)

# Plot SVGs with lower than average abundance in WM
feature = rownames(results_WM_both$low_genes)[seq_len(3)]
FeaturePlot(spe, feature, cluster_col = cluster_col, 
            coordinates = coordinates, cluster = 'WM', 
            legend_cluster = TRUE, annotation_cluster = TRUE, 
            linewidth = 0.6,title = TRUE)
