# Load required packages
library(data.table)    # For efficient data handling
library(tidyverse)    # For versatile data manipulation
library(ggplot2)       # For creating custom plots
library(ggrepel)       # For non-overlapping text labels
library(ggfortify)     # For statistical visualizations
library(ggprism)       # For publication-ready styling
library(pheatmap)      # For creating heatmaps
library(EnhancedVolcano) # For volcano plots
# Read differential expression results
deg_tbl <- fread("DEG_P1_lfc0.tsv")
dge_v <- readRDS("dge_v.rds")

# Set significance thresholds
FC_threshold <- 2
P_threshold <- 0.05

# Create color scheme for different gene categories
keyvals.colour_rd <- ifelse(
  (deg_tbl$logFC < -log2(FC_threshold) & deg_tbl$P.Value < P_threshold), 'blue',
  ifelse((deg_tbl$logFC > log2(FC_threshold) & deg_tbl$P.Value < P_threshold), 'red', 'grey30')
)

# Add informative labels showing number of genes in each category
names(keyvals.colour_rd)[keyvals.colour_rd == 'blue'] <- paste0(
  "Downregulated (n=", 
  sum(deg_tbl$logFC < -log2(FC_threshold) & deg_tbl$P.Value < P_threshold), 
  ")"
)

names(keyvals.colour_rd)[keyvals.colour_rd == 'red'] <- paste0(
  "Upregulated (n=", 
  sum(deg_tbl$logFC > log2(FC_threshold) & deg_tbl$P.Value < P_threshold), 
  ")"
)

names(keyvals.colour_rd)[keyvals.colour_rd == 'grey30'] <- 'Non-significant'
###
volcano_plot <- EnhancedVolcano(
  deg_tbl,
  lab = NA,
  x = 'logFC',
  y = 'adj.P.Val',
  title = '',
  subtitle = '',
  pCutoff = P_threshold,
  FCcutoff = log2(FC_threshold),
  gridlines.major = FALSE,
  gridlines.minor = FALSE,
  ylim = c(0, max(-log10(deg_tbl$adj.P.Val), na.rm=TRUE) + 0.5),
  colCustom = keyvals.colour_rd,
  axisLabSize = 20,
  labSize = 3,
  legendLabSize = 16,
  legendIconSize = 7,
  captionLabSize = 16,
  colAlpha = 1,
  pointSize = 0.3
)

### customize the volcano plot 
# Prepare data
deg_tbl$neg_log10_adj.P.Val <- -log10(deg_tbl$P.Value)
deg_tbl$significance <- ifelse(
  deg_tbl$P.Value < P_threshold & abs(deg_tbl$logFC) > log2(FC_threshold),
  ifelse(deg_tbl$logFC > 0, "Upregulated", "Downregulated"),
  "Non-Significant"
)

# Select top genes to label
genes_to_label <- deg_tbl[order(deg_tbl$P.Value)]$V1[1:15]

# Create enhanced volcano plot
custom_volcano <- ggplot(deg_tbl, aes(x = logFC, y = neg_log10_adj.P.Val)) +
  geom_point(aes(color = significance), alpha = 0.8, size = 1) +
  scale_color_manual(
    values = c("Upregulated" = "red", "Downregulated" = "blue", "Non-Significant" = "gray80")
  ) +
  geom_vline(xintercept = c(-log2(FC_threshold), log2(FC_threshold)), 
             linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = -log10(P_threshold), 
             linetype = "dashed", color = "gray50") +
  geom_text_repel(
    data = subset(deg_tbl, V1 %in% genes_to_label),
    aes(label = V1),
    size = 3,
    box.padding = 0.5,
    max.overlaps = Inf
  ) +
  labs(
    x = expression(log[2]~"Fold Change"),
    y = expression(-log[10]~"adjusted p-value"),
    color = "Differential Expression"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    text = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

