# Load required libraries
library(DESeq2)
library("pheatmap")
library("RColorBrewer")
library(ggplot2)
library(gridExtra)

# output directory
if (!dir.exists("plots"))
  dir.create("plots")

# Load the data ====
# Function to extract count data from the dataframe
extractCounts <- function(df, featureCols = NULL) {
  rownames(df) <- df$Geneid
  df <- df[, !(names(df)
               %in% unlist(c("Geneid", featureCols)))]
  return(df)
}

# Define the feature columns to extract and use later
featureCols <- c("Chrom", "Gene.Name", "Biotype", "Length")

# Load data
df <- read.csv("data/raw/ReverseStrandedCounts.csv")
cts <- extractCounts(df, featureCols = featureCols)
head(cts, 2)

# Create coldata
coldata <- data.frame(
  row.names = colnames(cts),
  condition = c(rep("treatment", 6), rep("control", 6)),
  timepoint = c(rep("3.5d", 3), rep("4.5h", 6), rep("3.5d", 3))
)
coldata$condition <-
  factor(coldata$condition, levels = c("control", "treatment"))
coldata$timepoint <-
  factor(coldata$timepoint, levels = c("4.5h", "3.5d"))
coldata <- coldata[order(coldata$timepoint, coldata$condition), ]
coldata

# Reorder counts matrix
cts <- cts[, rownames(coldata)]

# Define annotation colors for all plots generated thereafter
condition_colors <- c(control = "#A1C83E", treatment = "#F0978D")
timepoint_colors <- c(`4.5h` = "#63D7DE", `3.5d` = "#D79CF9")
annotation_colors = list(condition = condition_colors,
                         timepoint = timepoint_colors)

# Create the DESeq2 object ====
# Create DESeq dataset
dds <- DESeqDataSetFromMatrix(
  countData = cts,
  colData = coldata,
  design = ~ condition + timepoint + condition:timepoint
)
annotation_col <- as.data.frame(colData(dds))
dds <- DESeq(dds)

# Assign feature data to metadata columns
featureData <- data.frame(df[, featureCols])
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)

# Sample clustering and visualization ====
# Transform count data with rlog
rld <- rlog(dds, blind = TRUE)

# Heatmap of sample distances
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)

colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
p <- pheatmap(
  sampleDistMatrix,
  clustering_distance_rows = sampleDists,
  clustering_distance_cols = sampleDists,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  col = colors
)
ggsave(
  filename = "plots/sample_distance_heatmap.png",
  plot = p$gtable,
  device = "png",
  width = 8,
  height = 6
)

# Heatmap of genes
# Select the top 500 genes based on variance
gene_sd <- apply(counts(dds, normalized = TRUE), 1, var)
select <- order(gene_sd, decreasing = TRUE)[1:500]

p <- pheatmap(
  assay(rld)[select, ],
  show_rownames = FALSE,
  cluster_cols = TRUE,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  treeheight_row = 0
)
ggsave(
  filename = "plots/top_genes_heatmap.png",
  plot = p$gtable,
  device = "png",
  width = 6,
  height = 6
)

# PCA plot
pcaData <-
  plotPCA(rld,
          intgroup = c("condition", "timepoint"),
          returnData = TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))
p <-
  ggplot(pcaData, aes(PC1, PC2, color = condition, shape = timepoint)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  scale_color_manual(values = condition_colors) +
  scale_shape_manual(values = c(`4.5h` = 16, `3.5d` = 17))
ggsave(
  filename = "plots/pca_plot.png",
  plot = p,
  device = "png",
  width = 8,
  height = 6
)

# Differential expression analysis ====
# Extract condition effect specifically for each timepoint
dds$group <- factor(paste0(dds$timepoint, dds$condition))
design(dds) <- ~ group
dds <- DESeq(dds)

# 4.5 hour treatment vs. 4.5 hour control
res1 <-
  results(
    dds,
    contrast = c("group", "4.5htreatment", "4.5hcontrol"),
    alpha = 0.05,
    lfcThreshold = 1
  )

summary(res1)
png(
  filename = "plots/MA_plot_4.5h.png",
  width = 8,
  height = 6,
  units = "in",
  res = 72 * 2
)
plotMA(res1, ylim = c(-2, 2), alpha = 0.05)
dev.off()

# 3.5 day treatment vs. 3.5 day control
res2 <-
  results(
    dds,
    contrast = c("group", "3.5dtreatment", "3.5dcontrol"),
    alpha = 0.05,
    lfcThreshold = 1
  )

summary(res2)

png(
  filename = "plots/MA_plot_3.5d.png",
  width = 8,
  height = 6,
  units = "in",
  res = 72 * 2
)
p <- plotMA(res2, ylim = c(-2, 2), alpha = 0.05)
dev.off()

# Classify genes ====
# Significant genes for each timepoint
sig_genes_res1 <-
  rownames(res1)[which(res1$padj < 0.05 &
                         abs(res1$log2FoldChange) > 1)]
sig_genes_res2 <-
  rownames(res2)[which(res2$padj < 0.05 &
                         abs(res2$log2FoldChange) > 1)]

not_diff_expr <-
  setdiff(rownames(res1), union(sig_genes_res1, sig_genes_res2))

consistently_diff_expr <- intersect(sig_genes_res1, sig_genes_res2)

diff_expr_T1_only <- setdiff(sig_genes_res1, sig_genes_res2)

diff_expr_T2_only <- setdiff(sig_genes_res2, sig_genes_res1)

cat(
  "Number of genes consistently not differentially expressed:",
  length(not_diff_expr),
  "\n"
)
cat(
  "Number of genes consistently differentially expressed:",
  length(consistently_diff_expr),
  "\n"
)
cat("Number of genes differentially expressed only at T1:",
    length(diff_expr_T1_only),
    "\n")
cat("Number of genes differentially expressed only at T2:",
    length(diff_expr_T2_only),
    "\n")

p <- pheatmap(
  assay(rld[sig_genes_res1]),
  show_rownames = TRUE,
  show_colnames = FALSE,
  cluster_cols = FALSE,
  scale = "row",
  gaps_col = c(3, 6, 9),
  cutree_rows = 3,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  angle_col = 45
)
ggsave(
  filename = "plots/heatmap_diff_expr_T1.png",
  plot = p$gtable,
  device = "png",
  width = 8,
  height = 6
)


p <- pheatmap(
  assay(rld[sig_genes_res2, ]),
  show_rownames = FALSE,
  show_colnames = FALSE,
  cluster_cols = FALSE,
  scale = "row",
  gaps_col = c(3, 6, 9),
  cutree_rows = 2,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  angle_col = 45
)
ggsave(
  filename = "plots/heatmap_diff_expr_T2.png",
  plot = p$gtable,
  device = "png",
  width = 8,
  height = 6
)
