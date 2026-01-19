deg_tbl <- fread("DEG_P1_lfc0.tsv")
dge_v <- readRDS("dge_v.rds")
# violin
gene_to_plot <- "SGSM3"  # Replace with your gene of interest
dge_list <- sub("\\..*$", "", rownames(dge_v))
# Prepare data
library(biomaRt)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = dge_list,
  mart = mart
)
expr_mat <- as.data.frame(dge_v)
expr_mat$ensembl_gene_id <- sub("\\..*$", "", rownames(dge_v))
expr_mat <- expr_mat %>%
  left_join(gene_map, by = "ensembl_gene_id")
expr_mat <- expr_mat %>%
  filter(!is.na(hgnc_symbol), hgnc_symbol != "") %>%
  distinct(hgnc_symbol, .keep_all = TRUE)
rownames(expr_mat) <- expr_mat$hgnc_symbol
expr_mat$ensembl_gene_id <- NULL
expr_mat$hgnc_symbol <- NULL

expr_mat <- expr_mat + abs(min(expr_mat)) + 1
gene_vec <- expr_mat[gene_to_plot, ]   # expression values
gene_vec <- as.numeric(gene_vec)     
targets_ordered <- dge_v$targets[colnames(expr_mat), , drop = FALSE]
gene_data <- data.frame(
  Sample = colnames(expr_mat),
  Treatment = targets_ordered$condition,
  Expression = gene_vec
)

# Create violin plot
violin_plot <- ggplot(gene_data, aes(x = Treatment, y = Expression)) +
  geom_violin(aes(fill = Treatment), alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  stat_summary(fun = mean, geom = "point", size = 3, color = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  theme_minimal() +
  labs(title = paste(gene_to_plot, "Expression"), y = "log2 Expression") +
  theme(legend.position = "none")

# box plot
box_plot <- ggplot(gene_data, aes(x = Treatment, y = Expression)) +
  geom_boxplot(aes(fill = Treatment), alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  theme_minimal() +
  labs(title = paste(gene_to_plot, "Expression"), y = "log2 Expression") +
  theme(legend.position = "none")


# bar with SE
# Calculate means and standard errors
summary_data <- gene_data %>%
  group_by(Treatment) %>%
  summarise(
    mean = mean(Expression),
    se = sd(Expression)/sqrt(n())
  )

# Create bar plot
bar_plot <- ggplot(summary_data, aes(x = Treatment, y = mean, fill = Treatment)) +
  geom_bar(stat = "identity", alpha = 0.7) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 0.2) +
  theme_minimal() +
  labs(title = paste(gene_to_plot, "Expression"), y = "Mean log2 Expression") +
  theme(legend.position = "none")
