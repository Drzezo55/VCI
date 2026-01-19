library(data.table)
library(edgeR)
library(limma)
library(ggplot2)
library(ggrepel)
library(ggfortify)
library(stats)
library(sva)
#### create DGE 
dge <- DGEList(counts=df, samples = meta_c)
##Normalize
dge <- calcNormFactors(dge, method = "TMM")
dge_v <- voom(dge, plot=TRUE)
# save the dge_v object for later use
saveRDS(dge_v, "dge_v.rds")


# pca
shape_column <- "condition"
color_column <- "condition"
label <- TRUE
label_size <- 4
plot_save_name <- "PCA_Plot.png"

meta_table <- dge_v$targets
count_table_t <- as.data.frame(t(dge_v$E))
pca_prep <- prcomp(count_table_t, scale. = TRUE)

pca_plot <- autoplot(pca_prep, label, shape = shape_column, data = meta_table, colour = color_column) +
  geom_text_repel(aes(label = rownames(meta_table)), size = label_size) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank())

ggsave(plot_save_name, device = "png", units = "cm", width = 16, height = 14)
# 
comparison <- "VCI-Control"

design <- model.matrix(~ 0 + dge_v$targets[["condition"]])
colnames(design) <- gsub(".*]]", "", colnames(design))

contrast_matrix <- makeContrasts(contrasts = comparison, levels = design)

fit <- lmFit(dge_v, design)
fit_2 <- contrasts.fit(fit, contrast_matrix)
fit_2 <- eBayes(fit_2)

deg_tbl <- topTable(fit_2, coef = colnames(contrast_matrix), n = Inf, p.value=1, lfc=0, sort.by = "p")
fwrite(deg_tbl, "DEG_P1_lfc0.tsv", sep = "\t", row.names = T)

