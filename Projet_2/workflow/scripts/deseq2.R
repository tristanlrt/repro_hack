# =============================================================================
# Analyse différentielle par transcrit : DESeq2 1.26.0 (R 3.6)
# Objectif : reproduire la Fig. 2B de Hu et al. (2020), c'est-à-dire les
# transcrits différentiels (DETs) infecté vs contrôle à 11 et 33 dpi.
#
# Valeurs de l'article (section Résultats) :
#   11 dpi : 936 DETs (789 up, 147 down)
#   33 dpi : 2 081 DETs (1 908 up, 173 down)
#   communs aux deux temps : 618 up, 6 down
#
# Appelé par la règle deseq2 (workflow/regles/deseq2.smk) :
#   Rscript deseq2.R <samples.tsv> <sorties x5> <modele> <colonne_p>
#                    <seuil_p> <seuil_lfc> <quant.sf ...>
# =============================================================================

suppressPackageStartupMessages(library(DESeq2))
options(width = 160)

args <- commandArgs(trailingOnly = TRUE)
fichier_samples   <- args[1]
sortie_comptes    <- args[2]
sortie_d11        <- args[3]
sortie_d33        <- args[4]
sortie_comparaison <- args[5]
sortie_venn       <- args[6]
modele    <- args[7]               # "separe" ou "global"
colonne_p <- args[8]               # "pvalue" (brute) ou "padj" (ajustée BH)
seuil_p   <- as.numeric(args[9])   # 0.05 dans l'article
seuil_lfc <- as.numeric(args[10])  # 1.5 dans l'article (|log2FC| >= 1.5)
fichiers_quant <- args[11:length(args)]

stopifnot(modele %in% c("separe", "global"), colonne_p %in% c("pvalue", "padj"))

cat(R.version.string, "| DESeq2", as.character(packageVersion("DESeq2")), "\n")
cat("Modele principal :", modele, "| seuil :", colonne_p, "<", seuil_p,
    "et |log2FC| >=", seuil_lfc, "\n\n")


# -----------------------------------------------------------------------------
# 1. Matrice de comptes par transcrit (colonne NumReads des quant.sf)
# -----------------------------------------------------------------------------
# Hypothèse : on donne à DESeq2 les comptes estimés par Salmon (NumReads),
# arrondis à l'entier, sans agréger par gène (l'article teste des transcrits).
# Le nom de l'échantillon (SRR) est le nom du dossier qui contient quant.sf.
lire_quant <- function(fichier) {
  q <- read.delim(fichier, stringsAsFactors = FALSE)
  # gffread garde le préfixe "transcript:" du GFF3 Ensembl : on le retire
  setNames(q$NumReads, sub("^transcript:", "", q$Name))
}
quant <- lapply(fichiers_quant, lire_quant)
transcrits <- names(quant[[1]])
if (!all(vapply(quant, function(x) identical(names(x), transcrits), logical(1)))) {
  stop("Les quant.sf ne contiennent pas les mêmes transcrits dans le même ordre")
}
comptes <- round(do.call(cbind, quant))
storage.mode(comptes) <- "integer"
colnames(comptes) <- basename(dirname(fichiers_quant))
rownames(comptes) <- transcrits

write.table(data.frame(transcrit = rownames(comptes), comptes, check.names = FALSE),
            sortie_comptes, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Transcrits :", nrow(comptes), "| échantillons :", ncol(comptes), "\n")


# -----------------------------------------------------------------------------
# 2. Plan d'expérience (config/samples.tsv)
# -----------------------------------------------------------------------------
samples <- read.delim(fichier_samples, stringsAsFactors = FALSE)
rownames(samples) <- samples$sample
manquants <- setdiff(colnames(comptes), rownames(samples))
if (length(manquants) > 0) stop("Absents de samples.tsv : ", paste(manquants, collapse = ", "))
samples <- samples[colnames(comptes), ]
samples$infection <- factor(samples$infection, levels = c("controle", "infecte"))
samples$jour <- factor(samples$jour, levels = c("d11", "d33"))
samples$groupe <- factor(paste(samples$infection, samples$jour, sep = "_"),
                         levels = c("controle_d11", "infecte_d11", "controle_d33", "infecte_d33"))
print(table(samples$infection, samples$jour))


# -----------------------------------------------------------------------------
# 3. DESeq2 : deux façons de comparer infecté vs contrôle à chaque temps
# -----------------------------------------------------------------------------
# L'article ne donne ni le design ni les contrastes. On calcule les deux
# approches, et le choix principal (config deseq2: modele) sert aux tables
# et au diagramme de Venn. L'autre sert de comparaison dans le rapport.
#  - "separe" : un DESeq2 par temps, sur ses 6 échantillons (~ infection).
#               Correspond aux « two pair-wise comparisons » de l'article.
#  - "global" : un seul DESeq2 sur les 12 échantillons (~ groupe, 4 niveaux),
#               puis un contraste par temps. Équivalent au modèle avec
#               interaction ~ infection * jour pour ces deux contrastes.
lancer_deseq <- function(dds) {
  # La normalisation par défaut (médiane des ratios) exige au moins un
  # transcrit sans zéro. En mode test (chr19, 2 M de paires) ce n'est pas
  # garanti : on passe alors à "poscounts" et on le signale.
  if (any(rowSums(counts(dds) == 0) == 0)) {
    DESeq(dds, quiet = TRUE)
  } else {
    cat("ATTENTION : aucun transcrit sans zéro, normalisation 'poscounts'\n")
    DESeq(dds, sfType = "poscounts", quiet = TRUE)
  }
}

resultats <- list(separe = list(), global = list())
for (j in c("d11", "d33")) {
  garder <- samples$jour == j
  dds <- DESeqDataSetFromMatrix(comptes[, garder], samples[garder, ], design = ~ infection)
  resultats[["separe"]][[j]] <- results(lancer_deseq(dds),
                                        contrast = c("infection", "infecte", "controle"))
}
dds_global <- lancer_deseq(DESeqDataSetFromMatrix(comptes, samples, design = ~ groupe))
for (j in c("d11", "d33")) {
  resultats[["global"]][[j]] <- results(dds_global,
                                        contrast = c("groupe", paste0("infecte_", j), paste0("controle_", j)))
}


# -----------------------------------------------------------------------------
# 4. Sélection des DETs : p < seuil et |log2FC| >= seuil
# -----------------------------------------------------------------------------
# L'article est incohérent : les Méthodes disent « P-value < 0.05 » (brute),
# les Résultats « adjusted P-value < 0.05 ». Le tableau de comparaison
# donne les deux, la config choisit celle du Venn.
statut <- function(res, col) {
  p <- res[[col]]
  lfc <- res$log2FoldChange
  sig <- !is.na(p) & p < seuil_p & !is.na(lfc)
  ifelse(sig & lfc >= seuil_lfc, "up", ifelse(sig & lfc <= -seuil_lfc, "down", "non"))
}
liste_dets <- function(res, col, sens) rownames(res)[statut(res, col) == sens]


# -----------------------------------------------------------------------------
# 5. Tables complètes pour le modèle et le seuil choisis
# -----------------------------------------------------------------------------
ecrire_table <- function(res, fichier) {
  t <- data.frame(transcrit = rownames(res), as.data.frame(res),
                  DET = statut(res, colonne_p), stringsAsFactors = FALSE)
  t <- t[order(t$pvalue, na.last = TRUE), ]
  write.table(t, fichier, sep = "\t", quote = FALSE, row.names = FALSE)
}
ecrire_table(resultats[[modele]][["d11"]], sortie_d11)
ecrire_table(resultats[[modele]][["d33"]], sortie_d33)


# -----------------------------------------------------------------------------
# 6. Comparaison avec l'article (4 variantes : modèle x type de p-value)
# -----------------------------------------------------------------------------
compter <- function(r11, r33, col) {
  u11 <- liste_dets(r11, col, "up");   d11 <- liste_dets(r11, col, "down")
  u33 <- liste_dets(r33, col, "up");   d33 <- liste_dets(r33, col, "down")
  c(d11_total = length(u11) + length(d11), d11_up = length(u11), d11_down = length(d11),
    d33_total = length(u33) + length(d33), d33_up = length(u33), d33_down = length(d33),
    communs_up = length(intersect(u11, u33)), communs_down = length(intersect(d11, d33)))
}
lignes <- list(article = c(d11_total = 936, d11_up = 789, d11_down = 147,
                           d33_total = 2081, d33_up = 1908, d33_down = 173,
                           communs_up = 618, communs_down = 6))
for (m in c("separe", "global")) {
  for (col in c("pvalue", "padj")) {
    lignes[[paste(m, col, sep = "_")]] <- compter(resultats[[m]][["d11"]], resultats[[m]][["d33"]], col)
  }
}
comparaison <- data.frame(variante = names(lignes), do.call(rbind, lignes), row.names = NULL)
write.table(comparaison, sortie_comparaison, sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nComparaison avec l'article (variante principale :", paste(modele, colonne_p, sep = "_"), ")\n")
print(comparaison, row.names = FALSE)


# -----------------------------------------------------------------------------
# 7. Diagramme de Venn (Fig. 2B) : up et down, 11 dpi vs 33 dpi
# -----------------------------------------------------------------------------
# Dessiné en R de base (aucun paquet supplémentaire dans le conteneur).
# Sous chaque diagramme : les mêmes zones dans l'article, pour comparer.
venn2 <- function(a, b, titre, ref_article) {
  plot.new()
  plot.window(xlim = c(0, 10), ylim = c(0, 8), asp = 1)
  cercle <- function(x, col) {
    t <- seq(0, 2 * pi, length.out = 200)
    polygon(x + 2.6 * cos(t), 4 + 2.6 * sin(t), col = adjustcolor(col, 0.35), border = "grey30")
  }
  cercle(3.8, "#1b9e77")
  cercle(6.2, "#d95f02")
  text(c(2.4, 5, 7.6), 4, c(length(setdiff(a, b)), length(intersect(a, b)), length(setdiff(b, a))), cex = 1.3)
  text(c(3.2, 6.8), 7.2, c("11 dpi", "33 dpi"), font = 2)
  title(titre)
  mtext(paste0("Hu et al. 2020 : ", paste(ref_article, collapse = " | ")), side = 1, line = 0, cex = 0.8, col = "grey30")
}
res11 <- resultats[[modele]][["d11"]]
res33 <- resultats[[modele]][["d33"]]
pdf(sortie_venn, width = 9, height = 4.8)
par(mfrow = c(1, 2), mar = c(2, 1, 3, 1))
venn2(liste_dets(res11, colonne_p, "up"), liste_dets(res33, colonne_p, "up"),
      "Upregulated DETs", c(789 - 618, 618, 1908 - 618))
venn2(liste_dets(res11, colonne_p, "down"), liste_dets(res33, colonne_p, "down"),
      "Downregulated DETs", c(147 - 6, 6, 173 - 6))
invisible(dev.off())

cat("\n")
print(sessionInfo())
