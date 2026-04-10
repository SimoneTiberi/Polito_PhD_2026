Material for PhD course "Statistical Inference of Transcriptomics data" held at the Politecnico di Torino (Department of Mathematical Sciences) in March-April 2026.

Before the labs, install these packages:
- CRAN: ggplot2, ggfortify, ggforce, gplots, corrplot, UpSetR, purrr;
- Bioconductor: edgeR, EnhancedVolcano, clusterProfiler, DOSE, BANDITS, DEXSeq, fgsea, msigdbr, BiocParallel, ExperimentHub, scater, celldex, SingleR, scran, muscat, distinct, speckle, DESpace, SpatialExperiment, spatialLIBD, scuttle, BayesSpace.

For Bioconductor packages, first install BiocManager:
install.packages("BiocManager")
then install each package via:
BiocManager::install("edgeR")
BiocManager::install("EnhancedVolcano")
BiocManager::install("...")
