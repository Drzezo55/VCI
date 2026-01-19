# Prepare sample annotations
col_annot_df <- data.frame(
  Treatment = dge_v$targets$condition,
  row.names = rownames(dge_v$targets)
)
col_annot_df <- col_annot_df[order(col_annot_df$Treatment), , drop = FALSE]

# Select significant DEGs
deg_sig <- deg_tbl[abs(deg_tbl$logFC) > log2(FC_threshold) & 
                     deg_tbl$P.Value < P_threshold]$V1
expr_deg <- dge_v$E[deg_sig, rownames(col_annot_df), drop = FALSE]

# Create heatmap
heatmap_custom <- pheatmap(
  expr_deg,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  color = colorRampPalette(c("navy", "white", "red"))(100),
  annotation_col = col_annot_df,
  main = "Differential Expression Heatmap"
)
rownames(expr_deg) <- sub("\\..*$", "", rownames(expr_deg))

###
library(biomaRt)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = rownames(expr_deg),
  mart = mart
)
expr_df <- as.data.frame(expr_deg)
expr_df$ensembl_gene_id <- sub("\\..*$", "", rownames(expr_df))
expr_joined <- expr_df %>%
  left_join(gene_map, by = "ensembl_gene_id")
expr_joined <- expr_joined %>%
  filter(!is.na(hgnc_symbol), hgnc_symbol != "") %>%
  distinct(hgnc_symbol, .keep_all = TRUE)
expr_mat <- as.matrix(expr_joined[, colnames(expr_deg)])
rownames(expr_mat) <- expr_joined$hgnc_symbol
heatmap_custom <- pheatmap(
  expr_mat,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  color = colorRampPalette(c("navy", "white", "red"))(100),
  annotation_col = col_annot_df,
  main = "Differential Expression Heatmap"
)




save_pheatmap_pdf <- function(x, filename, width=8, height=8) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_pdf(heatmap_custom, "heatmap_custom.pdf", width=4, height=4)
