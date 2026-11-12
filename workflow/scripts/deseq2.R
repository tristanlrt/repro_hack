library(DESeq2)

args <- commandArgs(trailingOnly = TRUE)
fichier_counts <- args[1]
fichier_samples <- args[2]
fichier_genes <- args[3]
fichier_ribosome <- args[4]
fichier_trna <- args[5]
fichier_facteurs <- args[6]
sortie_table <- args[7]
sortie_resume <- args[8]
sortie_ma_complet <- args[9]
sortie_ma_traduction <- args[10]

counts <- read.table(fichier_counts, header = TRUE, sep = "\t", comment.char = "#", row.names = 1, check.names = FALSE)
counts <- counts[, 6:ncol(counts)]
colnames(counts) <- sub("\\.bam$", "", basename(colnames(counts)))
rownames(counts) <- sub("^gene-", "", rownames(counts))

samples <- read.table(fichier_samples, header = TRUE, sep = "\t", row.names = 1)
samples <- samples[colnames(counts), , drop = FALSE]
samples$condition <- factor(samples$condition, levels = c("control", "persister"))

dds <- DESeqDataSetFromMatrix(countData = as.matrix(counts), colData = samples, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds, alpha = 0.05)

genes <- read.delim(fichier_genes, header = FALSE, quote = "", stringsAsFactors = FALSE)
genes$id <- sub("^sao:", "", genes$V1)
genes$nom <- ifelse(grepl(";", genes$V4), sub(";.*", "", genes$V4), NA)

ribosome <- sub("^sao:", "", read.delim(fichier_ribosome, header = FALSE, stringsAsFactors = FALSE)$V2)
trna <- sub("^sao:", "", read.delim(fichier_trna, header = FALSE, stringsAsFactors = FALSE)$V2)
lignes_facteurs <- readLines(fichier_facteurs)
facteurs <- unlist(regmatches(lignes_facteurs, gregexpr("SAOUHSC_[0-9A-Za-z]+", lignes_facteurs)))
cds <- genes$id[genes$V2 == "CDS"]
traduction <- intersect(unique(c(ribosome, trna, facteurs)), cds)

table <- data.frame(gene = rownames(res), nom = genes$nom[match(rownames(res), genes$id)], as.data.frame(res), stringsAsFactors = FALSE)
table$traduction <- table$gene %in% traduction
write.table(table, sortie_table, sep = "\t", quote = FALSE, row.names = FALSE)

sig <- !is.na(res$padj) & res$padj < 0.05
trad <- table$traduction
resume <- c(
  paste("genes", nrow(res)),
  paste("genes testes", sum(!is.na(res$padj))),
  paste("DEG padj<0.05", sum(sig)),
  paste("up", sum(sig & res$log2FoldChange > 0)),
  paste("down", sum(sig & res$log2FoldChange < 0)),
  paste("genes traduction", sum(trad)),
  paste("genes traduction DEG", sum(trad & sig)),
  paste("genes traduction up", sum(trad & sig & res$log2FoldChange > 0)),
  paste("genes traduction down", sum(trad & sig & res$log2FoldChange < 0))
)
writeLines(resume, sortie_resume)

pdf(sortie_ma_complet, width = 7, height = 6)
plotMA(res, alpha = 0.05, ylim = c(-4, 4), main = "MA-plot persisters vs control")
dev.off()

t <- table[trad & !is.na(table$baseMean) & table$baseMean > 0, ]
t_sig <- !is.na(t$padj) & t$padj < 0.05
t_facteur <- t$gene %in% facteurs
t_trna <- t$gene %in% trna
pdf(sortie_ma_traduction, width = 7, height = 6)
plot(log2(t$baseMean), t$log2FoldChange, pch = 20, col = ifelse(t_sig, "red", "grey50"),
     xlim = c(0, 20), ylim = c(-6, 5), xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n",
     xlab = "log2 base Mean", ylab = "log2 Fold Change", main = "Genes de la traduction")
axis(1, at = seq(0, 20, 2))
axis(2, at = -6:5, las = 1)
abline(h = 0, lty = 2)
points(log2(t$baseMean[t_trna]), t$log2FoldChange[t_trna], cex = 1.3)
text(log2(t$baseMean[t_facteur]), t$log2FoldChange[t_facteur], t$nom[t_facteur], pos = 3, cex = 0.7)
legend("bottomleft", c("padj < 0.05", "non significatif", "AA-tRNA synthetases"), pch = c(20, 20, 1), col = c("red", "grey50", "black"), bty = "n")
dev.off()
