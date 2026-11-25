# Projet 2 : cerveau de souris infectée par *Toxoplasma gondii*

Reproduction partielle de Hu et al. (2020), *Front. Microbiol.* 11:570903.
Workflow Snakemake, où chaque étape tourne dans un conteneur Apptainer (Damona/Zenodo, Biocontainers).
Objectif : la Fig. 2B, diagramme de Venn des transcrits différentiels (DETs) à 11 et 33 dpi.

```
SRA (PRJNA483261) ──> Trimmomatic 0.39 ──┐
                                         ├──> salmon quant ──> quant.sf ──> DESeq2 1.26.0 ──> Venn (Fig. 2B)
Ensembl 100 (GRCm38) ─> gffread ─> salmon index (sans leurres)
```

## Lancer le test (chromosome 19)

```bash
cd Projet_2
./install.sh          # une seule fois : Apptainer + Snakemake 9.27.1 (conda)
./run.sh -n           # dry-run : affiche les tâches sans rien exécuter (43 tâches)
./run.sh              # exécution
```

La config par défaut ([config/config.yaml](config/config.yaml)) est en **mode test** :

| paramètre    | test      | analyse complète |
|--------------|-----------|------------------|
| `chromosome` | `"19"`    | `""` (toplevel)  |
| `max_spots`  | `2000000` | `0` (toutes les lectures) |

Pour tester un seul échantillon :
`./run.sh results/salmon/chr19_sub2000000/SRR7657878/quant.sf`

Résultats :
- `results/trimmed/<lect>/` : lectures nettoyées (`_1P` / `_2P`)
- `results/index/salmon_<ref>/` : index Salmon
- `results/salmon/<ref>_<lect>/<SRR>/quant.sf` : quantification par transcrit
- `results/deseq2/<ref>_<lect>/` : résultats DESeq2 (voir ci-dessous)
- `logs/` : logs Trimmomatic, Salmon et DESeq2. Le taux d'assignation est dans
  `results/salmon/.../aux_info/meta_info.json` (`percent_mapped`). Sur chr19
  seul, il est faible (~3-5 %), ce qui est attendu.

## Analyse différentielle (DESeq2)

Règle `deseq2` ([workflow/regles/deseq2.smk](workflow/regles/deseq2.smk)), script
[workflow/scripts/deseq2.R](workflow/scripts/deseq2.R). Elle attend les 12 `quant.sf`, puis :

1. matrice de comptes par transcrit : colonne `NumReads` des `quant.sf`, arrondie (DESeq2 exige des entiers) ;
2. plan d'expérience lu dans `config/samples.tsv` (infecté / contrôle × 11 / 33 dpi, 3 souris par groupe) ;
3. DESeq2 : infectées vs contrôles à chaque temps ;
4. DETs : p < 0,05 et |log2FC| ≥ 1,5, séparés en up / down ;
5. diagramme de Venn 11 dpi vs 33 dpi, pour les up et pour les down.

Sorties dans `results/deseq2/<ref>_<lect>/` :

| Fichier | Contenu |
|---|---|
| `comparaison_article.tsv` | nombres de DETs pour 4 variantes (modèle × p-value), à côté de ceux de l'article |
| `venn_dets.pdf` | Fig. 2B reproduite, avec les chiffres de l'article sous chaque diagramme |
| `resultats_d11.tsv`, `resultats_d33.tsv` | résultats DESeq2 complets + colonne `DET` (up / down / non) |
| `comptes_transcrits.tsv` | matrice de comptes |

Paramètres dans la section `deseq2:` de la config : `modele` (`separe` ou `global`),
`colonne_p` (`pvalue` ou `padj`), `seuil_p`, `seuil_log2fc`.

Chiffres de l'article à retrouver : 936 DETs à 11 dpi (789 up, 147 down), 2 081 à 33 dpi
(1 908 up, 173 down), dont 618 up et 6 down communs.

## Plan d'expérience (reconstruit depuis le run selector SRA)

Voir [config/samples.tsv](config/samples.tsv). Les codes T111…C333 sont ceux du
Table 1 de l'article, déduits du champ `LibraryName` de SRA
(ex. « test group 1 on 11 day post-infection » → T111).

## Hypothèses et lacunes de reproductibilité

| Point | Ce que dit l'article | Choix fait | Justification |
|---|---|---|---|
| Type de lectures | rien | paired-end 2×150 | `LibraryLayout = PAIRED`, `avgLength = 300` sur SRA |
| Séquenceur | HiSeq 2500 | — | SRA base de données indique **HiSeq X Ten** (incohérence) |
| Lectures déposées | — | — | le nombre de bases sur SRA est égal aux *clean bases* du Table 1 : les FASTQ déposés semblent **déjà nettoyés** par Novogene, donc Trimmomatic devrait retirer peu de lectures |
| Paramètres Trimmomatic | « adaptor and low-quality reads removed » | `ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:2:True LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36` | paramètres par défaut du manuel ; NEBNext Ultra utilise des adaptateurs compatibles TruSeq |
| Transcriptome | « GRCm38.p6, Ensembl 100 » | transcrits extraits du génome + GFF3 avec gffread 0.12.7 | consigne du projet (construire le cDNA à partir du génome + annotation) |
| Index Salmon | rien | k = 31, **sans leurres** (`decoys: false`) | l'article ne parle pas de leurres ; avec le génome entier en leurre, l'index demande beaucoup plus de mémoire et de temps. Les leurres restent activables dans la config |
| Type de librairie | rien | `-l A` (auto-détection) | NEBNext Ultra (non directionnel) : détecté `IU` en test |
| Options quant | « bias-aware quantification » | `--seqBias --gcBias --validateMappings` | l'article mentionne explicitement la correction de biais |
| Version DESeq2 | « DEseq2 package (v.3.6.1) » | DESeq2 1.26.0 sous R 3.6.1 | 3.6.1 n'est pas une version de DESeq2 (toutes en 1.x) : probablement la version de R. Bioconductor 3.10, conçu pour R 3.6.1, contient DESeq2 1.26.0 |
| Entrée de DESeq2 | rien | `NumReads` de Salmon arrondis, par transcrit | la Fig. 2B compte des transcrits, pas des gènes |
| Design DESeq2 | ni design ni contrastes | `separe` par défaut (un DESeq2 par temps, 6 souris) ; `global` (12 souris, `~ groupe`) calculé aussi | rien ne permet de trancher : `comparaison_article.tsv` compare les deux |
| Seuil de p-value | Méthodes : P < 0.05 ; Résultats : *adjusted* P < 0.05 | p brute par défaut, `padj` calculé aussi | incohérence de l'article : les deux sont comparés |
| Seuil de log2FC | Méthodes : ≥ 1.5 ; Résultats : > 1.5 | ≥ 1,5 | Méthodes ; l'écart ne joue que sur une égalité exacte |

## Conteneurs utilisés

| Outil | Version | Image |
|---|---|---|
| sra-tools | 3.1.0 | Damona (Zenodo 10894351) |
| wget | 1.24.5 | Damona (Zenodo 14415567) |
| Trimmomatic | 0.39 | Damona (Zenodo 22145522) |
| Salmon | 1.2.0 | Damona (Zenodo 22144681) |
| gffread | 0.12.7 | Damona (Zenodo 8037039) |
| DESeq2 | 1.26.0 (R 3.6.1) | Biocontainers `bioconductor-deseq2:1.26.0--r36he1b5a44_0` (référencée par Damona) ; recette équivalente : [containers/deseq2-1.26/Dockerfile](containers/deseq2-1.26/Dockerfile) |

## Environnements testés et pièges

Mode test (chr19) exécuté de bout en bout (43 tâches, ~1 h 50) sur une VM IFB
(Ubuntu, 4 cœurs, 15 Go de RAM), avec Apptainer 1.5.5 et Snakemake **9.27.2**
(environnement `snakemake` déjà présent sur la VM, donc pas la 9.27.1 de `install.sh`).

- **Lancer depuis le vrai chemin du dossier.** Sur les VM IFB, `~/data` est un lien vers
  `/ifb/data`. Lancé depuis `~/data/...`, `fastq-dump` échoue dans le conteneur
  (`directory unauthorized`) : utiliser `cd /ifb/data/mydatalocal/.../Projet_2`.
- **Scripts exécutables.** Si `./install.sh` répond `Permission denied` : `chmod +x install.sh run.sh`.
- **Mode complet.** Prévoir de la place : les FASTQ complets des 12 échantillons et leurs
  versions nettoyées occupent plusieurs dizaines de Go.
