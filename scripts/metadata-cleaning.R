df <- read.csv("data/GSE303449_raw.csv")
View(df)
meta <- read.delim("data/GSE303449_series_matrix.txt",
  header = FALSE,
  quote = "\"",
  sep = "\t",
  stringsAsFactors = FALSE
)
meta <- t(meta)
View(meta)
meta <- as.data.frame(meta)
library(DataEditR)
meta <- data_edit(meta)
colnames(meta) <- meta[1,]
meta <- meta[-1,]
# Remove the prefix from every row
meta$sex <- sub("^Sex: ", "", meta$sex)
meta$condition <- sub("^condition: ", "", meta$condition)
meta$race <- sub("^race: ", "", meta$race)
meta$pmi <- sub("^pmi: ", "", meta$pmi)
meta$age <- sub("^age: ", "", meta$age)
meta$pmi <- as.numeric(meta$pmi)
meta$age <- as.numeric(meta$age)
colnames(df[,-1]) %in% meta$id
# Fix column names in your expression matrix
colnames(df)[-1] <- gsub("\\.", "-", colnames(df)[-1])
colnames(df[,-1]) %in% meta$id
colnames(df[,-1]) %in% meta$id
##
df$genes <- df$ID
df$ID <- NULL
rownames(df) <- df$genes
df$genes <- NULL
meta$condition <- as.factor(meta$condition)
rownames(meta) <- meta$id

#creat metadata
meta_c <- meta[, c("id", "condition")]
#### Filter 
perc_keep <- 0.8
gene_keep <- rowSums(df > 0) >= ceiling(perc_keep * ncol(df))
df <- df[gene_keep, ]
