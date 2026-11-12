library(DESeq2)

args <- commandArgs(trailingOnly = TRUE)
fichier_counts <- args[1]
fichier_samples <- args[2]
fichier_genes <- args[3]
fichier_noms <- args[4]
fichier_ribosome <- args[5]
fichier_trna <- args[6]
fichier_facteurs <- args[7]
sortie_table <- args[8]
sortie_resume <- args[9]
sortie_ma_complet <- args[10]
sortie_ma_traduction <- args[11]

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

noms <- read.delim(fichier_noms, quote = "", stringsAsFactors = FALSE)
noms$nom <- ifelse(noms[, 2] == "-", noms[, 3], noms[, 2])

ribosome <- sub("^sao:", "", read.delim(fichier_ribosome, header = FALSE, stringsAsFactors = FALSE)$V2)
trna <- sub("^sao:", "", read.delim(fichier_trna, header = FALSE, stringsAsFactors = FALSE)$V2)
lignes_facteurs <- readLines(fichier_facteurs)
facteurs <- unlist(regmatches(lignes_facteurs, gregexpr("SAOUHSC_[0-9A-Za-z]+", lignes_facteurs)))
cds <- genes$id[genes$V2 == "CDS"]
traduction <- intersect(unique(c(ribosome, trna, facteurs)), cds)

table <- data.frame(gene = rownames(res), nom = noms$nom[match(rownames(res), noms[, 1])], as.data.frame(res), stringsAsFactors = FALSE)
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

ok <- res$baseMean > 0 & !is.na(res$log2FoldChange)
m <- res$baseMean[ok]
lfc <- res$log2FoldChange[ok]
s <- sig[ok]
haut <- lfc > 4
bas <- lfc < -4
lfc <- pmax(pmin(lfc, 4), -4)
couleur <- ifelse(s, rgb(1, 0, 0, 0.6), rgb(0, 0, 0, 0.6))
forme <- ifelse(haut, 17, ifelse(bas, 25, 16))

pdf(sortie_ma_complet, width = 7, height = 6)
par(mar = c(5, 5, 2, 2))
plot(m, lfc, log = "x", type = "n", ylim = c(-4.2, 4.2), xaxt = "n", yaxt = "n", bty = "n",
     xlab = "Mean of normalized counts", ylab = expression(log[2] ~ "fold change"))
u <- par("usr")
rect(10^u[1], u[3], 10^u[2], u[4], col = "grey92", border = NA)
abline(v = 10^seq(-1, 5, 1), col = "grey97", lwd = 0.5)
abline(h = seq(-3, 3, 2), col = "grey97", lwd = 0.5)
abline(v = 10^seq(0, 4, 2), h = seq(-4, 4, 2), col = "white", lwd = 1.2)
abline(h = 0, lty = 2, col = "grey40")
points(m, lfc, pch = forme, col = couleur, bg = couleur, cex = ifelse(haut | bas, 0.6, 0.4))
axis(1, at = 10^c(0, 2, 4), labels = expression(10^0, 10^2, 10^4), col = NA, col.ticks = "grey40")
axis(2, at = seq(-4, 4, 2), las = 1, col = NA, col.ticks = "grey40")
dev.off()

t <- table[trad & !is.na(table$baseMean) & table$baseMean > 0, ]
t_sig <- !is.na(t$padj) & t$padj < 0.05
t_trna <- t$gene %in% trna
x <- log2(t$baseMean)
y <- t$log2FoldChange

etiquettes <- data.frame(
  nom = c("frr", "infA", "tsf", "infC", "infB", "pth"),
  dx = c(-1.5, -1, 2, 3, 3, -4),
  dy = c(2, 2.5, 1.5, 0, -0.5, -0.5),
  stringsAsFactors = FALSE
)

pdf(sortie_ma_traduction, width = 6, height = 6)
par(mar = c(5, 5, 2, 2))
plot(x, y, type = "n", xlim = c(0, 20), ylim = c(-6, 5), xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n",
     xlab = expression(Log[2] ~ "base Mean"), ylab = expression(Log[2] ~ "Fold Change"))
axis(1, at = seq(0, 20, 2))
axis(2, at = -6:5, las = 1)
abline(h = 0, lty = 2, lwd = 2)
points(x, y, pch = 16, col = ifelse(t_sig, "red", "grey50"))
points(x[t_trna], y[t_trna], pch = 1, cex = 1.2, lwd = 1.5)
for (i in seq_len(nrow(etiquettes))) {
  j <- which(t$nom == etiquettes$nom[i])
  if (length(j) == 1) {
    xt <- x[j] + etiquettes$dx[i]
    yt <- y[j] + etiquettes$dy[i]
    segments(x[j], y[j], xt, yt, lwd = 1.5)
    text(xt, yt, etiquettes$nom[i], font = 3, pos = ifelse(etiquettes$dx[i] < 0, 2, 4), offset = 0.2)
  }
}
legend("bottomleft", c("Significant", "Non-significant"), pch = 16, col = c("red", "grey50"), bty = "n", cex = 0.8)
legend("bottomright", "AA-tRNA synthetases", pch = 1, pt.lwd = 1.5, text.font = 3, bty = "n", cex = 0.8)
box(lwd = 2)
dev.off()
