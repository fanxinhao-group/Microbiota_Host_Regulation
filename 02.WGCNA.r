library(WGCNA)
library(reshape2)
library(stringr)
library(matrixStats)
library(tidyverse)
library(corrplot)
library(tidygraph)

library(dplyr)
library(tidyr)
library(igraph)
library(uwot)
library(ggraph)
library(ggplot2)
library(ggrepel)
library(data.table)

options(stringsAsFactors = FALSE)
# Enable multithreading.
enableWGCNAThreads()

# The recommended network types are "signed" and "signed hybrid".

type = "signed hybrid"
# Correlation calculation.
# The recommended correlation method is biweight midcorrelation (bicor).
# corType: pearson or bicor

corType = "bicor"

setwd("\\WGCNA_path\\")

##all tissues
# TPM file: the first column contains gene IDs, followed by sample columns.
all_tpm <- read.csv("gene_tpm_matrix.csv",header = T,row.names = 1)

# Trait file: contains sample groups or phenotypes; SampleIDs must match the TPM column names.
trait <- read.csv("trait.csv", header = TRUE,row.names = 1)
trait$gas_ratio <- trait$gastrocnemius_weight / trait$body_weight * 100
trait$sol_ratio <- trait$soleus_weight / trait$body_weight * 100
trait$TA_ratio <- trait$TA_weight / trait$body_weight * 100
trait$quad_ratio <- trait$quadriceps_weight / trait$body_weight * 100
trait <- trait[,c("GF","SPF","CON","gas_ratio","TA_ratio","quad_ratio")]

## Filter lowly expressed genes.
keep <- rowSums(tpm > 1) >= ncol(tpm) * 0.5
tpm_filt <- tpm[keep, ]

## Apply log2 transformation.
tpm_log <- log2(tpm_filt + 1)
## Select highly variable genes using the median absolute deviation (MAD).
gene_mad <- rowMads(as.matrix(tpm_log))
cutoff <- quantile(gene_mad, 0.5)

tpm_var <- tpm_log[gene_mad > cutoff, ]
## Transpose the matrix so that samples are rows and genes are columns.
dataExpr <- as.data.frame(t(tpm_var))

nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)

## Check data quality.
gsg <- goodSamplesGenes(dataExpr, verbose = 3)

## Identify potential outlier samples.
pdf("Sample_clustering_to_detect_outliers.pdf",10,6)
sampleTree = hclust(dist(dataExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
dev.off()

pdf("Soft_Threshold_power.pdf",10,6)
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(dataExpr, powerVector=powers, 
                        networkType=type, verbose=5)
par(mfrow = c(1,2))
cex1 = 0.85

# The x-axis represents the soft-thresholding power; the y-axis shows the scale-free topology fit. Higher values indicate a better fit.
# A higher fit index indicates a network that more closely follows scale-free topology.
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# Selection criterion: scale-free topology fit index (R-squared) = 0.85.
abline(h=0.85,col="red")
# Relationship between the soft-thresholding power and mean connectivity.
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers,
     cex=cex1, col="red")
dev.off()

power = sft$powerEstimate
power

## One-step network construction and module detection.

net <- blockwiseModules(
  dataExpr,
  power = power,
  networkType = "signed hybrid",
  TOMType = "signed",
  corType = "bicor",
  maxPOutliers = 0.1,
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = FALSE,
  verbose = 3
)
table(net$colors)
save(net,file = "WGCNA_net.RData")

## Genes assigned to the grey module are not assigned to any detected co-expression module.
# Convert labels to colors for plotting
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
modules <- cbind(as.data.frame(moduleLabels),as.data.frame(moduleColors))

write.csv(modules,"modules.csv",quote=F)
# Visualize the dendrogram with module colors displayed underneath.
# If the results are unsatisfactory, recutBlockwiseTrees() can be used to adjust the module boundaries while saving computation time.
pdf("cluster_dendrogram.pdf",10,6)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

## Module–trait associations.
# Calculate module eigengenes (MEs).
MEs <- moduleEigengenes(dataExpr, moduleColors)$eigengenes
MEs <- orderMEs(MEs)
MEs
write.table(MEs, "module_eigengenes.txt", sep="\t", quote=FALSE)
# Plot correlations among modules based on gene expression patterns.
# The marDendro and marHeatmap arguments set the bottom, left, top, and right margins, respectively.
pdf("plotEigengeneNetworks.pdf",10,8)
plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)
dev.off()

## Calculate correlations.
moduleTraitCor <- cor(MEs, trait, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(dataExpr))
## Generate the heatmap.
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep="")
save(moduleTraitCor,moduleTraitPvalue,textMatrix,file = "Module-trait_ratio_relationships.RData")

text_matrix <- matrix(paste(signif(moduleTraitCor, 2), "\n(",
                            signif(moduleTraitPvalue, 1), ")", sep=""), nrow = nrow(moduleTraitCor))
rownames(text_matrix) <- rownames(moduleTraitCor)
colnames(text_matrix) <- colnames(moduleTraitCor)

# Extract module colors.
module_names_clean <- sub("^ME", "", rownames(moduleTraitCor))

# Create module annotations.
row_anno <- rowAnnotation(
  Module = module_names_clean,
  col = list(Module = setNames(module_names_clean, module_names_clean)),
  show_legend = c(Module = FALSE),
  annotation_name_gp = gpar(fontsize = 12, fontface = "bold")
)
# Create the heatmap.
ht <- Heatmap(
  moduleTraitCor,
  name = "cor",
  column_title = "Module−trait relationships",
  column_title_gp = gpar(fontsize = 20,col = "#2C3E50"),
  column_title_side = "top",
  left_annotation = row_anno,
  row_names_side = "left",
  show_row_names = T,
  show_column_names = T,
  cluster_rows = F,
  cluster_columns = F,
  cell_fun = function(j, i, x, y, width, height, fill) {
    label <- text_matrix[i, j] 
    grid.text(label, x, y, gp = gpar(fontsize = 12, col = "black"))},
  row_names_gp = gpar(fontsize = 15),
  column_names_gp = gpar(fontsize = 15),
  column_names_rot = 30,
  column_names_max_height = unit(4, "cm"),
  heatmap_legend_param = list(
    color_bar = "continuous",     
    legend_height = unit(8, "cm"), 
    legend_width = unit(1, "cm")  
  )
)
pdf("Module-trait_relationships_ratio.pdf",10,8)
draw(ht)
dev.off()
