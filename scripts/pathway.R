remove.packages(c("ggplot2", "data.table"))
install.packages('Rcpp', dependencies = TRUE)
install.packages('ggplot2', dependencies = TRUE)
install.packages('data.table', dependencies = TRUE)
library(data.table)
library(clusterProfiler)
# library(org.Hs.eg.db)
library(org.Hs.eg.db)
# library(org.Rn.eg.db)
library(ggplot2)
library(DOSE)
library(enrichplot)
deg <- fread("DEG_P1_lfc0.tsv")
p_threshold <- 2
fc_threshold <- 2
# rank the DEGs by the fold change
deg_order_fc <- deg[order(-logFC)] # rank the genes by logFC in descending order
logfc <- deg_order_fc$logFC # get logFC
names(logfc) <- sub("\\..*$", "", deg_order_fc$V1)

# GSEA
enrich_go_gsea <- gseGO(
  geneList  = logfc,
  OrgDb     = org.Hs.eg.db,
  ont       = "ALL",       # BP / CC / MF / ALL
  keyType   = "ENSEMBL",    
  pvalueCutoff = 0.9,
  verbose   = FALSE
)
# save the enrichment table
enrich_go_gsea_df <- enrich_go_gsea@result

fwrite(enrich_go_gsea_df, "enrich_go_gsea_df.tsv", sep = "\t")
dotplot_enrich_go_gsea <- dotplot(enrich_go_gsea, showCategory = 10, orderBy="GeneRatio")

ggsave("dotplot_enrich_go_gsea.png", dotplot_enrich_go_gsea, device = "png", units = "cm", width = 16, height = 18)

##


# specify the GO term you want to plot (check your "enrich_go_gsea_df" above)
GOID <- "GO:0042060"
gsea_plot <- gseaplot2(enrich_go_gsea, geneSetID = GOID, title = subset(enrich_go_gsea_df, ID == GOID)$Description)

ggsave(paste0("gsea_plot_", gsub(":", "", GOID), ".png"), gsea_plot, device = "png", units = "cm", width = 18, height = 16)

### ORA 
# make a gene list
deg_up <- deg[logFC > log2(fc_threshold) & adj.P.Val < p_threshold]$V1
deg_dn <- deg[logFC < -log2(fc_threshold) & adj.P.Val < p_threshold]$V1
deg_up <- sub("\\..*$", "", deg_up)
deg_dn <- sub("\\..*$", "", deg_dn)
enrich_go_fet_up <- enrichGO(gene = deg_up, OrgDb=org.Hs.eg.db, keyType="ENSEMBL", ont="ALL", pvalueCutoff=0.9, pAdjustMethod="BH", qvalueCutoff=0.9)
enrich_go_fet_up_df <- enrich_go_fet_up@result
fwrite(enrich_go_fet_up_df, "enrich_go_fet_up_df.tsv", sep = "\t")

# show the top 10 terms
dotplot_enrich_go_fet_up <- dotplot(enrich_go_fet_up, showCategory = 10, orderBy="GeneRatio")

ggsave("dotplot_enrich_go_fet_up.png", dotplot_enrich_go_fet_up, device = "png", units = "cm", width = 16, height = 18)

##
summary(deg$logFC)
summary(deg$adj.P.Val)
table(deg$logFC > log2(fc_threshold))
table(deg$adj.P.Val < p_threshold)
##




