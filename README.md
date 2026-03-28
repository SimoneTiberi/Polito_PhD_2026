Material for PhD course "Statistical Inference of Transcriptomics data" held at the Politecnico di Torino in 2026.

Before the labs, install these packages:
- CRAN: ggplot2, ggfortify, gplots, corrplot, UpSetR, purrr;
- Bioconductor: edgeR, EnhancedVolcano, clusterProfiler, DOSE, BANDITS, DEXSeq, fgsea, msigdbr, BiocParallel, ExperimentHub, scater, celldex, SingleR, scran, muscat, distinct, speckle.

For Bioconductor packages, first install BiocManager:
install.packages("BiocManager")
then install each package via:
BiocManager::install("edgeR")
BiocManager::install("EnhancedVolcano")
BiocManager::install("...")
