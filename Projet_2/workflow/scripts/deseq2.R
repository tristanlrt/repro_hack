# =============================================================================
# Analyse différentielle par transcrit : DESeq2 1.26.0 (R 3.6)
# Objectif : reproduire la Fig. 2B de Hu et al. (2020), c'est-à-dire les
# transcrits différentiels (DETs) infecté vs contrôle à 11 et 33 dpi.
#
# Valeurs de l'article (section Résultats et Fig. 2B) :
#   11 dpi : 936 DETs (789 up, 147 down)
#   33 dpi : 2 081 DETs (1 908 up, 173 down)
#   communs aux deux temps : 618 up, 6 down ; sens inversé : 2 + 2
#
# Appelé par la règle deseq2 (workflow/regles/deseq2.smk) :
#   Rscript deseq2.R <samples.tsv> <sorties x6> <modele> <colonne_p>
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
sortie_pca        <- args[7]
modele    <- args[8]               # "separe" ou "global"
colonne_p <- args[9]               # "pvalue" (brute) ou "padj" (ajustée BH)
seuil_p   <- as.numeric(args[10])  # 0.05 dans l'article
seuil_lfc <- as.numeric(args[11])  # 1.5 dans l'article (|log2FC| >= 1.5)
fichiers_quant <- args[12:length(args)]

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
#               L'article ne précise pas le design : les deux sont calculés.
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
    communs_up = length(intersect(u11, u33)), communs_down = length(intersect(d11, d33)),
    up11_down33 = length(intersect(u11, d33)), down11_up33 = length(intersect(d11, u33)))
}
lignes <- list(article = c(d11_total = 936, d11_up = 789, d11_down = 147,
                           d33_total = 2081, d33_up = 1908, d33_down = 173,
                           communs_up = 618, communs_down = 6,
                           up11_down33 = 2, down11_up33 = 2))
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
# 7. Diagramme de Venn (Fig. 2B) : 4 ensembles, comme dans l'article
# -----------------------------------------------------------------------------
# Les 4 listes (up et down, à 11 et 33 dpi) sur un seul Venn à 4 ellipses :
# il montre aussi les transcrits qui changent de sens entre les deux temps.
# À gauche nos chiffres, à droite ceux de la Fig. 2B, aux mêmes positions.
# Dessiné en R de base (aucun paquet supplémentaire dans le conteneur).
noms <- c("Acute_up", "Acute_down", "Chronic_up", "Chronic_down")
ellipses <- data.frame(x = c(0.35, 0.47, 0.53, 0.65), y = c(0.45, 0.55, 0.55, 0.45),
                       a = 0.36, b = 0.19, angle = c(-50, -50, 50, 50))
dans_ellipse <- function(px, py, e) {
  t <- -e$angle * pi / 180
  dx <- px - e$x; dy <- py - e$y
  ((dx * cos(t) - dy * sin(t)) / e$a)^2 + ((dx * sin(t) + dy * cos(t)) / e$b)^2 <= 1
}
# Position de chaque zone : point de la zone le plus proche de son centre
grille <- expand.grid(x = seq(0, 1, length.out = 301), y = seq(0, 1, length.out = 301))
code <- rep(0, nrow(grille))
for (k in 1:4) code <- code + (2^(k - 1)) * dans_ellipse(grille$x, grille$y, ellipses[k, ])
positions <- t(sapply(1:15, function(z) {
  g <- grille[code == z, ]
  cx <- mean(g$x); cy <- mean(g$y)
  unlist(g[which.min((g$x - cx)^2 + (g$y - cy)^2), ])
}))
# Code de zone = somme des ensembles : Acute_up 1, Acute_down 2, Chronic_up 4, Chronic_down 8
compter_zones <- function(listes) {
  tous <- unique(unlist(listes))
  z <- rep(0, length(tous))
  for (k in 1:4) z <- z + (2^(k - 1)) * (tous %in% listes[[k]])
  tabulate(z, nbins = 15)
}
article_zones <- rep(0, 15)
article_zones[c(1, 2, 4, 8, 5, 10, 6, 9)] <- c(169, 139, 1288, 165, 618, 6, 2, 2)
venn4 <- function(zones, titre) {
  plot.new()
  plot.window(xlim = c(-0.05, 1.05), ylim = c(0, 1.05), asp = 1)
  couleurs <- c("#1f78b4", "#e6c229", "#e31a1c", "#33a02c")
  for (k in 1:4) {
    e <- ellipses[k, ]
    t <- seq(0, 2 * pi, length.out = 200)
    px <- e$a * cos(t); py <- e$b * sin(t); r <- e$angle * pi / 180
    polygon(e$x + px * cos(r) - py * sin(r), e$y + px * sin(r) + py * cos(r),
            col = adjustcolor(couleurs[k], 0.25), border = "grey30")
  }
  text(positions[, "x"], positions[, "y"], zones, cex = 0.9)
  text(c(0.08, 0.33, 0.67, 0.92), c(0.93, 1.0, 1.0, 0.93), noms, cex = 0.75, font = 2)
  title(titre, cex.main = 1)
}
res11 <- resultats[[modele]][["d11"]]
res33 <- resultats[[modele]][["d33"]]
nos_zones <- compter_zones(list(liste_dets(res11, colonne_p, "up"), liste_dets(res11, colonne_p, "down"),
                                liste_dets(res33, colonne_p, "up"), liste_dets(res33, colonne_p, "down")))
pdf(sortie_venn, width = 11, height = 5.5)
par(mfrow = c(1, 2), mar = c(1, 1, 3, 1))
venn4(nos_zones, paste0("Ce workflow (", modele, ", ", colonne_p, " < ", seuil_p, ")"))
venn4(article_zones, "Hu et al. 2020, Fig. 2B")
invisible(dev.off())

# -----------------------------------------------------------------------------
# 8. PCA des 12 échantillons (contrôle qualité, absente de l'article)
# -----------------------------------------------------------------------------
# L'article affirme une « clear separation between infected and control
# transcriptomes » sans la montrer : la PCA permet de le vérifier.
# Comptes stabilisés par vst() en mode aveugle (sans utiliser les groupes),
# puis plotPCA() de DESeq2 sur les 500 transcrits les plus variables.
# En mode test (chr19), trop peu de transcrits pour vst() : on passe à
# varianceStabilizingTransformation(), équivalente mais plus lente.
vsd <- tryCatch(vst(dds_global, blind = TRUE),
                error = function(e) varianceStabilizingTransformation(dds_global, blind = TRUE))
pca <- plotPCA(vsd, intgroup = c("infection", "jour"), returnData = TRUE)
part <- round(100 * attr(pca, "percentVar"))
pca$code <- samples[rownames(pca), "code"]
couleur <- ifelse(pca$infection == "infecte", "#d7301f", "#2c7fb8")
forme <- ifelse(pca$jour == "d11", 16, 17)
# Textes de la figure sans accents : le conteneur tourne en locale "C".
pdf(sortie_pca, width = 8, height = 6, useDingbats = FALSE)
par(mar = c(5, 5, 3, 8))
marge_x <- diff(range(pca$PC1)) * 0.15
marge_y <- diff(range(pca$PC2)) * 0.15
plot(pca$PC1, pca$PC2, col = couleur, pch = forme, cex = 1.8,
     xlim = range(pca$PC1) + c(-marge_x, marge_x), ylim = range(pca$PC2) + c(-marge_y, marge_y),
     xlab = paste0("PC1 (", part[1], "% of variance)"),
     ylab = paste0("PC2 (", part[2], "% of variance)"),
     main = "PCA of the 12 samples (vst, 500 most variable transcripts)", cex.main = 0.95)
text(pca$PC1, pca$PC2, pca$code, pos = 3, cex = 0.75, col = "grey30")
# Légende hors du cadre, pour ne jamais cacher de points
legend(par("usr")[2], par("usr")[4], xpd = NA, bty = "n", cex = 0.85,
       legend = c("Infected", "Control", "11 dpi", "33 dpi"),
       col = c("#d7301f", "#2c7fb8", "grey30", "grey30"), pch = c(15, 15, 16, 17))
invisible(dev.off())
cat("\nCoordonnées de la PCA :\n")
write.table(data.frame(code = pca$code, sample = rownames(pca), PC1 = round(pca$PC1, 2), PC2 = round(pca$PC2, 2)),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n")
print(sessionInfo())
