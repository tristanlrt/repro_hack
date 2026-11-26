# repro_hack
projet reprohackathon

# Projet 1 : persisters intracellulaires de *Staphylococcus aureus* — Peyrusson et al. (2020)

Reproduction partielle de Peyrusson F. et al. (2020), *Intracellular Staphylococcus aureus persisters upon antibiotic exposure*, Nature Communications 11:2200.

Workflow Snakemake dans lequel chaque étape s'exécute dans un conteneur Apptainer. Les versions des outils sont reproduites autant que possible à l'identique de celles de l'article.

## En bref

Notre pipeline retrouve **1 482 gènes différentiellement exprimés (711 up, 771 down)**, contre **1 477 (710 up, 767 down)** dans l'article, soit un écart de seulement 0,3 %.

- **Article :** Peyrusson F. et al. (2020), *Intracellular Staphylococcus aureus persisters upon antibiotic exposure*, Nature Communications 11:2200.
- **À reproduire :** la Fig. 3c (MA-plots des gènes liés à la traduction) et la Supplementary Fig. 3 (MA-plot de tous les gènes).
- **Comparaison :** persisters intracellulaires (bactéries présentes dans des macrophages J774 exposées à l'oxacilline) contre contrôles (bactéries extracellulaires mélangées à un lysat de cellules J774).
- **Données :** 6 échantillons RNA-seq provenant de GEO (GSE139659), correspondant à 3 réplicats persisters et 3 réplicats contrôles.
- **Dépôt :** [github.com/tristanlrt/repro_hack](https://github.com/tristanlrt/repro_hack), branche `main`.

Le pipeline comprend 8 étapes, de la récupération des données SRA jusqu'à la production des MA-plots et des résultats DESeq2.

## Lancer le pipeline

Le pipeline est conçu pour être exécuté sur une machine Linux x86_64, par exemple une VM IFB sous Ubuntu 24.04.

Trois commandes suffisent :

```bash
git clone https://github.com/tristanlrt/repro_hack.git
cd repro_hack

./install.sh       # Installe Apptainer, Miniforge et Snakemake si absents
./run.sh           # Lance l'ensemble du pipeline
```

L'installation prépare :
- Apptainer 1.5.5 ;
- Miniforge si conda n'est pas disponible ;
- un environnement conda `snakemake` avec Snakemake 9.27.1 et pandas.

Le script `install.sh` peut être relancé sans risque : il n'installe que les éléments manquants.

### Exécution sur une VM

Le pipeline prend plusieurs heures. Il est donc recommandé de le lancer dans une session `tmux`, qui continue à fonctionner même après la déconnexion SSH.

```bash
tmux new -s pipeline

./run.sh 2>&1 | tee pipeline.log
```

Pour se détacher de la session sans arrêter le calcul : `Ctrl+b`, puis `d`.

Pour retrouver la session :

```bash
tmux attach -t pipeline
```

### Temps d'exécution

| Machine | Durée depuis zéro |
|---|---:|
| VM 4 vCPU / 16 Go (IFB ifb.m4.xlarge), mesurée | ≈ 4 h 45 |
| VM 16 vCPU / 16 Go (IFB ifb.tr.2xlarge), estimée | ≈ 2 h |

Les résultats sont générés dans `results/figures/` (MA-plots PDF) et `results/deseq2/` (tableaux et résumé).

## Structure du dépôt

Le dépôt ne contient que le code et les fichiers de configuration. Les données et résultats sont téléchargés ou générés à l'exécution et ne sont pas versionnés (`.gitignore`).

| Fichier / dossier | Rôle |
|---|---|
| `install.sh` | Prépare l'environnement d'exécution : Apptainer, Miniforge et Snakemake. |
| `run.sh` | Active l'environnement conda et lance Snakemake avec Apptainer. |
| `config/config.yaml` | Paramètres du pipeline : threads, URLs des références, ressources KEGG, AureoWiki et images des outils. |
| `config/samples.tsv` | Métadonnées des 6 échantillons et conditions expérimentales. |
| `containers/*/Dockerfile` | Recettes de construction des 4 images personnalisées : cutadapt, bowtie, subread et DESeq2. |
| `workflow/Snakefile` | Point d'entrée du workflow : chargement de la configuration, inclusion des règles et définition des résultats finaux. |
| `workflow/regles/telechargement.smk` | Téléchargement des FASTQ, du génome, de l'annotation et des ressources KEGG/AureoWiki. |
| `workflow/regles/index.smk` | Construction de l'index Bowtie du génome. |
| `workflow/regles/trimming.smk` | Nettoyage des lectures avec cutadapt. |
| `workflow/regles/mapping.smk` | Alignement avec Bowtie, tri et indexation des BAM avec samtools. |
| `workflow/regles/comptage.smk` | Comptage des lectures par gène avec featureCounts. |
| `workflow/regles/deseq2.smk` | Exécution de l'analyse différentielle. |
| `workflow/scripts/deseq2.R` | Analyse DESeq2, annotation des gènes, calcul des statistiques et génération des trois MA-plots. |

### Dossiers créés à l'exécution

Ces dossiers ne sont pas versionnés :

- `data/raw/` : les 6 fichiers FASTQ téléchargés (environ 9 Go).
- `data/ref/` : génome de référence, annotation et ressources KEGG/AureoWiki.
- `results/index/` : index Bowtie.
- `results/trimmed/` : lectures nettoyées.
- `results/bam/` : fichiers BAM alignés, triés et indexés.
- `results/counts/` : matrices de comptage.
- `results/deseq2/` : résultats de l'analyse différentielle.
- `results/figures/` : MA-plots au format PDF.
- `.snakemake/` : état interne de Snakemake et cache des images Apptainer.

## Le pipeline étape par étape

Huit étapes permettent de passer des 6 runs SRA aux résultats d'expression différentielle. Depuis zéro, le workflow comporte **31 jobs Snakemake**.

### Schéma général

```text
                    DONNÉES RNA-seq
                    6 runs SRA
                         |
                         v
                   FASTQ × 6
                         |
                         v
                  Trimming × 6
                    cutadapt
                         |
                         v
                  Mapping × 6  <----- Génome de référence
                    Bowtie             CP000253.1
                         |                  |
                         v                  v
                    BAM × 6           Index Bowtie
                         |
                         v
                     Comptage  <----- Annotation GFF
                   featureCounts
                         |
                         v
                    counts.txt
                         |
                         v
                      DESeq2  <----- KEGG + AureoWiki
                         |
                         v
                 Résultats différentiels
                         |
              +----------+----------+
              |          |          |
              v          v          v
           MA-plot    MA-plot    MA-plot
           complet   traduction  traduction
                       stricte     élargie
```

Le génome, l'annotation et les ressources d'annotation fonctionnelle sont téléchargés indépendamment des lectures, puis rejoignent le workflow au moment du mapping, du comptage et de l'analyse différentielle.

### Détail des étapes

| Étape (règle) | Outil / commande | Ce qu'il faut savoir |
|---|---|---|
| 1. `telecharger_fastq` | `fasterq-dump --threads 4 --progress` puis gzip | Téléchargement des 6 runs SRA. Lectures single-end de 65 nt, entre 27 et 52 millions de reads par échantillon. |
| 2. `telecharger_genome` et `telecharger_annotation` | `wget` | Téléchargement du génome de référence et de son annotation GFF. |
| 3. `indexer_genome` | `bowtie-build` | Construction de l'index Bowtie. Le génome ne fait que 2,8 Mo : cette étape prend quelques secondes. |
| 4. `trimmer` | `cutadapt -q 20 -m 25` | Filtrage des lectures selon leur qualité et leur longueur minimale. Moins de 0,002 % des reads sont écartés. |
| 5. `mapper` | `bowtie -p 4 -S` puis `samtools sort` et `samtools index` | Alignement des lectures sur le génome de référence et génération des fichiers BAM. |
| 6. `compter` | `featureCounts` | Comptage des lectures par gène à partir des 6 BAM et de l'annotation GFF. |
| 7. `telecharger_kegg` et `telecharger_aureowiki` | API REST KEGG et requêtes AureoWiki | Récupération des listes de gènes liés à la traduction et des noms des gènes. |
| 8. `deseq2` | `Rscript workflow/scripts/deseq2.R` | Analyse différentielle, annotation des résultats et génération des MA-plots. |

### Ce que fait le script DESeq2

Le script `workflow/scripts/deseq2.R` réalise les opérations suivantes :

1. Lit `counts.txt`, conserve les 6 colonnes de comptage et retire le préfixe `gene-` des identifiants.
2. Lit `config/samples.tsv` et définit `control` comme condition de référence. Un log2FC positif correspond donc à une expression plus élevée chez les persisters.
3. Construit un objet `DESeqDataSetFromMatrix` avec le design `~ condition`.
4. Exécute `DESeq()`, puis `results(alpha = 0.05)` : test de Wald, correction de Benjamini-Hochberg et filtrage indépendant.
5. Construit deux listes de gènes liés à la traduction, uniquement à partir des gènes codants (CDS) :
   - **Liste stricte (92 gènes)** : ribosome, aminoacyl-ARNt synthétases et facteurs de traduction.
   - **Liste élargie (182 gènes)** : mêmes catégories, auxquelles s'ajoutent la biogenèse du ribosome et celle des ARNt.
6. Génère les tableaux de résultats, le résumé chiffré et les trois MA-plots.

## Les échantillons

Les 6 échantillons proviennent du jeu de données GEO GSE139659.

| Identifiant SRA | Condition | Réplicat |
|---|---|---|
| SRR10379721 | Persister | 1 |
| SRR10379722 | Persister | 2 |
| SRR10379723 | Persister | 3 |
| SRR10379724 | Control | 1 |
| SRR10379725 | Control | 2 |
| SRR10379726 | Control | 3 |

Les métadonnées sont définies dans `config/samples.tsv`.

Les trois premiers échantillons correspondent aux bactéries persisters intracellulaires, tandis que les trois derniers constituent les contrôles extracellulaires.

## Les conteneurs utilisés

Chaque étape du pipeline est exécutée dans un conteneur afin de garantir la reproductibilité des versions logicielles.

Damona ne proposait pas les versions exactes de quatre outils utilisés dans l'article. Nous avons donc construit nos propres images pour cutadapt, Bowtie, Subread et DESeq2.

| Outil | Version | Image / origine |
|---|---|---|
| sra-tools (`fasterq-dump`) | 3.1.0 | Damona, Zenodo |
| wget | 1.24.5 | Damona, Zenodo |
| cutadapt | 1.11 | `docker://tristanlrt/cutadapt:1.11`, construite depuis Bioconda |
| Bowtie | 0.12.7 | `docker://tristanlrt/bowtie:0.12.7`, binaire officiel SourceForge |
| samtools | 1.9 | Même image que Bowtie, depuis Bioconda |
| Subread (`featureCounts`) | 1.4.6-p3 | `docker://tristanlrt/subread:1.4.6-p3`, binaire officiel SourceForge |
| R + DESeq2 | R 3.4.1 + DESeq2 1.16.1 | `docker://tristanlrt/deseq2:1.16` |

### Construction des images personnalisées

- **cutadapt :** installation de `cutadapt=1.11` depuis Bioconda via micromamba.
- **Bowtie :** utilisation du binaire officiel 0.12.7 depuis SourceForge, associé à samtools 1.9.
- **Subread :** utilisation du binaire officiel précompilé 1.4.6-p3 depuis SourceForge, cette version n'étant pas disponible sur Bioconda.
- **DESeq2 :** utilisation de `r-base=3.4.1`, `bioconductor-deseq2=1.16.1` et `libgfortran=3`.

La dépendance `libgfortran=3` est indispensable : sans elle, R ne démarre pas correctement et renvoie une erreur concernant `libgfortran.so.3`.

Les images ont été construites sur une VM Linux x86_64 puis publiées sur Docker Hub. Cette étape est importante : un Mac Apple Silicon produirait des images ARM incompatibles avec les machines Linux x86_64 utilisées pour l'exécution.

Snakemake télécharge ensuite automatiquement les images et les convertit au format Apptainer.

## Résultats

Les résultats obtenus sont très proches de ceux de l'article. Le nombre total de gènes différentiellement exprimés ne diffère que de 5 gènes.

### Comparaison avec l'article

| Indicateur | Notre reproduction | Article |
|---|---:|---:|
| Gènes annotés | 2 967 | 2 967 |
| Gènes différentiellement exprimés (`padj < 0,05`) | 1 482 | 1 477 |
| Up chez les persisters | 711 | 710 |
| Down chez les persisters | 771 | 767 |

Sur les 2 967 gènes annotés, 2 752 ont été testés par DESeq2. Les autres ne possèdent aucun read ou ont été écartés par le filtrage indépendant.

DESeq2 normalise les différences de profondeur de séquençage entre les échantillons.

### Qualité des données par échantillon

| Échantillon | Condition | Reads (millions) | Écartés par cutadapt | Alignement | Reads assignés à un gène (millions) |
|---|---|---:|---:|---:|---:|
| SRR10379721 | Persister 1 | 51,9 | 658 | 63,4 % | 24,3 |
| SRR10379722 | Persister 2 | 44,2 | 463 | 61,5 % | 20,4 |
| SRR10379723 | Persister 3 | 39,7 | 287 | 61,8 % | 18,3 |
| SRR10379724 | Control 1 | 47,3 | 442 | 67,3 % | 22,3 |
| SRR10379725 | Control 2 | 46,2 | 512 | 77,2 % | 26,2 |
| SRR10379726 | Control 3 | 27,0 | 303 | 36,1 % | 6,4 |

Le contrôle 3 (SRR10379726) présente un taux d'alignement plus faible (36,1 %) que les autres échantillons (61 à 77 %). La cause de cette différence n'a pas été étudiée.

### Gènes liés à la traduction — Fig. 3c

Deux listes de gènes ont été utilisées pour reproduire la figure de l'article.

| Liste | Nombre de gènes | Différentiellement exprimés | Up | Down |
|---|---:|---:|---:|---:|
| Stricte (`ma_plot_traduction.pdf`) | 92 | 43 | 24 | 19 |
| Élargie (`ma_plot_traduction_large.pdf`) | 182 | 93 | 54 | 39 |
| Article (estimation à l'œil) | ≈ 190 | ≈ 100 | — | — |

Les deux listes permettent de retrouver la conclusion biologique de l'article :

- Les aminoacyl-ARNt synthétases sont majoritairement réprimées chez les persisters.
- Plusieurs facteurs de traduction, notamment `frr`, `infA`, `infB`, `infC` et `tsf`, sont activés.

La liste élargie se rapproche davantage de la figure publiée, avec environ 180 points et 93 gènes différentiellement exprimés.

### Fichiers produits

| Fichier | Contenu |
|---|---|
| `results/figures/ma_plot_complet.pdf` | Reproduction de la Supplementary Fig. 3 : MA-plot de tous les gènes. |
| `results/figures/ma_plot_traduction.pdf` | Reproduction de la Fig. 3c avec la liste stricte de 92 gènes. |
| `results/figures/ma_plot_traduction_large.pdf` | Reproduction de la Fig. 3c avec la liste élargie de 182 gènes. |
| `results/deseq2/resultats_deseq2.tsv` | Résultats DESeq2 complets, avec les noms des gènes et leur appartenance aux deux listes. |
| `results/deseq2/resume.txt` | Résumé des principaux résultats et statistiques. |

Les MA-plots reproduisent les axes, les couleurs, les légendes et les étiquettes de l'article. Pour les figures de traduction, les six mêmes gènes sont étiquetés : `frr`, `infA`, `tsf`, `infC`, `infB` et `pth`.

## Choix et hypothèses méthodologiques

L'article ne précise pas tous les paramètres nécessaires à la reproduction. Plusieurs choix ont donc dû être effectués et justifiés.

| Question | Notre choix | Justification |
|---|---|---|
| Retrait des adaptateurs | Commande du sujet sans option `-a` | L'article indique que les séquences adaptatrices ont été retirées, mais ne donne pas leur séquence. Les tests montrent que les reads contiennent quasiment aucun adaptateur. |
| Seuil de qualité | `-q 20` | Valeur proposée dans le sujet, correspondant à une qualité Phred 20 (environ 1 erreur pour 100 bases). |
| Orientation de la librairie | `-s 1` | Vérifiée expérimentalement sur 200 000 reads : 47 509 assignés avec `-s 1`, contre 47 837 avec `-s 0` et seulement 578 avec `-s 2`. |
| Référence DESeq2 | `control` | Un log2FC positif correspond à une expression plus élevée chez les persisters, conformément au sens de l'article. |
| Seuil DESeq2 | `padj < 0,05` | L'article utilise la p-value ajustée et le filtrage indépendant. Le paramètre `alpha = 0.05` aligne le filtrage sur ce seuil. |
| Noms des gènes | AureoWiki, souche NCTC8325 | Le GFF NCBI fournit principalement des identifiants de locus. AureoWiki permet de récupérer les noms des gènes. |
| Liste des gènes de traduction | Deux listes : stricte (92) et élargie (182) | L'article ne détaille pas précisément les catégories KEGG retenues. Les deux listes permettent de comparer les résultats. |
| ARNt et ARNr | Exclus des listes de traduction | La figure de l'article montre des protéines. Les ARN ajoutaient de nombreux points peu exprimés qui ne semblaient pas présents dans la figure. |
| Étiquettes de la Fig. 3c | Les six mêmes gènes que l'article | `frr`, `infA`, `tsf`, `infC`, `infB`, `pth`, afin de faciliter la comparaison visuelle. |
| Version samtools | 1.9 | Version disponible sur Bioconda et compatible avec la commande de mapping du sujet. |
| Style des figures | R de base | Le conteneur ne contient que les versions de R et DESeq2 de l'article, sans ggplot2. |
| Outils d'exécution | Snakemake 9.27.1 et Apptainer 1.5.5 | Apptainer fonctionne sans droits root et permet l'exécution des conteneurs sur la VM. |

### Construction des listes de traduction

L'article mentionne les gènes liés à la traduction, classés selon leur appartenance à la KEGG Orthology, mais ne précise pas toutes les catégories retenues.

Nous avons donc construit deux listes :

**Liste stricte (92 gènes)**

- Ribosome (KEGG 03010).
- Aminoacyl-ARNt synthétases (KEGG 00970).
- Facteurs de traduction (KEGG 03012).

Elle correspond à la catégorie « Translation » de KEGG.

**Liste élargie (182 gènes)**

Elle reprend les trois catégories précédentes et ajoute :

- Biogenèse du ribosome (KEGG 03009).
- Biogenèse des ARNt (KEGG 03016).

Cette deuxième liste est celle qui se rapproche le plus de la figure publiée : environ 180 points et 93 gènes différentiellement exprimés, contre environ 190 et 100 dans l'article.

L'écart restant pourrait être lié à l'évolution de KEGG depuis 2018, année de la référence utilisée par les auteurs. Cette hypothèse n'a toutefois pas pu être vérifiée.

## Écarts avec les commandes du sujet

Les commandes proposées dans le sujet ont été reprises, à l'exception de trois adaptations nécessaires.

| Commande du sujet | Adaptation dans notre pipeline | Justification |
|---|---|---|
| `fasterq-dump <SRAID>` puis gzip `*.fastq` | Téléchargement et compression échantillon par échantillon | Les échantillons sont traités en parallèle. Une compression globale pourrait interférer avec les fichiers des autres échantillons encore en cours d'écriture. |
| `featureCounts --extraAttributes Name ...` | Option supprimée | Cette option n'existe pas dans Subread 1.4.6-p3. Les noms des gènes sont récupérés séparément depuis AureoWiki. |
| Fichiers directement dans le dossier courant | Organisation dans `results/`, `data/ref/` et `data/raw/` | Organisation des fichiers par étape pour faciliter la reproductibilité et la maintenance. Les options biologiques ne changent pas. |

## Fonctionnement sur la VM

Lorsque `./run.sh` est exécuté, Snakemake construit automatiquement le graphe des dépendances à partir des fichiers finaux définis dans `rule all`.

### Déroulement

1. **Activation :** `run.sh` active l'environnement conda contenant Snakemake.
2. **Lecture :** Snakemake lit `workflow/Snakefile`, `config/config.yaml` et `config/samples.tsv`, puis inclut les différentes règles.
3. **Construction du graphe :** Snakemake remonte les dépendances depuis les résultats finaux pour déterminer les étapes à exécuter.
4. **Gestion des conteneurs :** les images sont téléchargées depuis Zenodo ou Docker Hub, converties au format Apptainer et mises en cache dans `.snakemake/`.
5. **Exécution :** chaque job est lancé dans son conteneur avec le dossier du projet monté dans celui-ci.
6. **Parallélisation :** Snakemake répartit les jobs selon le nombre de cœurs disponibles. Les étapes mono-thread, comme cutadapt 1.11, peuvent traiter plusieurs échantillons en parallèle.
7. **Fin :** si toutes les étapes réussissent, le journal indique `N of N steps (100%) done`.

Snakemake ne relance que les étapes dont les sorties sont absentes, plus anciennes que leurs entrées ou affectées par une modification du code ou des paramètres.

En cas d'échec, il est possible de relancer `./run.sh` : Snakemake reprend à partir des étapes qui restent à effectuer.

### Commandes utiles

```bash
# Afficher les tâches qui seraient exécutées sans les lancer
./run.sh -n

# Recalculer DESeq2 et les figures
./run.sh --forcerun deseq2

# Suivre le journal d'exécution
tail -f pipeline.log

# Consulter le résumé des résultats
cat results/deseq2/resume.txt
```

Les résultats restent stockés sur le disque de la VM. Pour récupérer les figures depuis son ordinateur :

```bash
scp "<utilisateur>@<ip-de-la-vm>:repro_hack/results/figures/*.pdf" .
```

**Attention : supprimer la VM entraîne la suppression des résultats qui n'ont pas été récupérés.**

## Temps de calcul et ressources

Le pipeline complet a été exécuté sur une VM IFB de 4 vCPU et 16 Go de RAM.

| Étape | Durée mesurée | Principal facteur limitant |
|---|---:|---|
| Téléchargement des 6 FASTQ et compression gzip | 2 h 02 | Réseau NCBI et compression mono-thread. |
| Trimming, mapping et comptage | 2 h 42 | cutadapt 1.11 mono-thread et mappings successifs. |
| Génome, annotation, index, KEGG, AureoWiki, DESeq2 et figures | Moins d'une minute au total | — |
| **Total** | **≈ 4 h 45** | **Processeur et réseau.** |

Ressources utilisées :

- Mémoire : environ 1 Go.
- Espace disque : environ 25 Go.
- FASTQ : jusqu'à 9,4 Go décompressés par fichier.

Sur une VM de 16 vCPU, la durée est estimée à environ 2 heures. Le gain reste toutefois limité par les étapes mono-thread, notamment gzip et cutadapt.

## Limites et points ouverts

Plusieurs éléments restent à discuter ou pourraient être approfondis :

- **Liste des gènes de traduction :** l'article ne fournit pas la liste exacte. Deux listes ont été construites, mais le choix définitif reste à discuter.
- **Adaptateurs :** les adaptateurs n'ont pas été retirés. Les tests indiquent qu'ils sont quasiment absents des données, mais ce choix reste une différence par rapport à la description de l'article.
- **Échantillon SRR10379726 :** le contrôle 3 présente un taux d'alignement de 36,1 %, nettement inférieur aux autres échantillons. La cause n'a pas été identifiée.
- **Services externes :** le pipeline dépend de NCBI, KEGG, AureoWiki, Zenodo et Docker Hub. Les ressources KEGG et AureoWiki peuvent évoluer dans le temps, ce qui peut légèrement modifier les résultats lors d'une nouvelle exécution.
- **Écart avec l'article :** cinq gènes différentiellement exprimés séparent notre résultat de celui de la publication. Les causes exactes n'ont pas été vérifiées. Les adaptateurs, les versions de samtools et de l'annotation ou certains paramètres non publiés pourraient contribuer à cet écart.
- **Reproductibilité des figures :** les deux listes de traduction permettent de comparer les résultats, mais la liste exacte utilisée par les auteurs n'est pas connue.

## Conclusion

Le pipeline permet de reproduire les principales analyses transcriptomiques de Peyrusson et al. (2020), depuis les données RNA-seq brutes jusqu'à l'analyse différentielle et aux MA-plots.

La reproduction retrouve 1 482 gènes différentiellement exprimés, contre 1 477 dans l'article, soit une concordance très élevée.

L'utilisation de Snakemake et Apptainer permet d'automatiser les étapes, de contrôler les versions des logiciels et de faciliter la réexécution du workflow sur une machine Linux.

Les principaux écarts concernent la définition des gènes liés à la traduction et certains paramètres méthodologiques insuffisamment détaillés dans l'article. Ces choix sont documentés afin de rendre la démarche transparente et reproductible.