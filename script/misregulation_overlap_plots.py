# TODO take from [[id:59330f09-5bfa-4be9-9a5e-c52f2a81870e][Revision plots for miRNA target prediction paper]]
library(eulerr)

# hacky
df.col.readnames <- c("Geneid", "Genename",
                      "Drosha_tpm", "Drosha_log2FoldChange", "Drosha_padj",
                      "Dicer_tpm", "Dicer_log2FoldChange", "Dicer_padj",
                      "Ago12_tpm", "Ago12_log2FoldChange", "Ago12_padj",
                      "Ago1_tpm", "Ago1_log2FoldChange", "Ago1_padj",
                      "Ago2_tpm", "Ago2_log2FoldChange", "Ago2_padj",
                      "WT_2i_tpm", "WT_2i_log2FoldChange", "WT_2i_padj",
                      "WT_tpm", "WT_log2FoldChange", "WT_padj")

df.col.usenames <- c("Geneid", "Genename",
                      "Drosha_tpm", "Drosha_log2FoldChange", "Drosha_padj",
                      "Dicer_tpm", "Dicer_log2FoldChange", "Dicer_padj",
                      "Ago12_tpm", "Ago12_log2FoldChange", "Ago12_padj",
                      "WT_tpm", "WT_log2FoldChange", "WT_padj")

df <- read.table(snakemake@input[[1]], sep=",", skip=2, col.names=df.col.readnames)
df <- df[df.col.usenames]  # ignore ago single mutants

mutants <- snakemake@params[["mutants"]]

# only select genes with minimal expression (say 0.5 tpm)
mcols <- sapply(mutants, function(m) paste(m, "tpm", sep="_"))
max_values <- apply(df[, mcols], 1, FUN=max)
df <- df[max_values > 0.5, ]
df <- df[complete.cases(df), ]  # drop nans

# generate a matrix with true false values for our columns:
log2fccols <- sapply(mutants, function(m) paste(m, "log2FoldChange", sep="_"))
padjcols <- sapply(mutants, function(m) paste(m, "padj", sep="_"))
if (snakemake@wildcards[["direction"]] == "up") {
  mat <- (df[, log2fccols] > 0.5) & (df[, padjcols] < 0.05)
} else {
  mat <- (df[, log2fccols] < -0.5) & (df[, padjcols] < 0.05)
}
colnames(mat) <- mutants

fit <- euler(mat)

svg(snakemake@output[[1]])
plot(fit,
     alpha = 0.7,
     fills = list(fill = c("#009292", "#FFB6DB", "#B66DFF", "#6DB6FF")),
     quantities = list(cex=2.0),
     labels=list(cex=2.5))  # cex is the font size
dev.off()
