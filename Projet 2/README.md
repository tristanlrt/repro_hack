# Projet 2 : cerveau de souris infectée par *Toxoplasma gondii*

Reproduction partielle de Hu et al. (2020), *Front. Microbiol.* 11:570903.
Workflow Snakemake, où chaque étape tourne dans un conteneur Apptainer (Damona/Zenodo).

```
SRA (PRJNA483261) ──> Trimmomatic 0.39 ──┐
                                         ├──> salmon quant ──> quant.sf (par transcrit)
Ensembl 100 (GRCm38) ─> gffread ─> salmon index (decoys)
```

## Lancer le test (chromosome 19)

```bash
cd "Projet 2"
./install.sh          # une seule fois : Apptainer + Snakemake 9.27.1 (conda)
./run.sh -n           # dry-run : affiche les tâches sans rien exécuter
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
- `logs/` : logs Trimmomatic et Salmon. Le taux d'assignation est dans
  `results/salmon/.../aux_info/meta_info.json` (`percent_mapped`). Sur chr19
  seul, il est faible (~3-5 %), ce qui est attendu.

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
| Index Salmon | rien | k = 31, génome en leurres (decoys) | valeur par défaut pour lectures longues ; les decoys sont recommandés par les auteurs de Salmon |
| Type de librairie | rien | `-l A` (auto-détection) | NEBNext Ultra (non directionnel) : on s'attend à `IU` |
| Options quant | « bias-aware quantification » | `--seqBias --gcBias --validateMappings` | l'article mentionne explicitement la correction de biais |
| Version DESeq2 | « v3.6.1 » | à décider | ce n'est pas une version de DESeq2 (probablement R 3.6.1) |

## Conteneurs utilisés

| Outil | Version | Image |
|---|---|---|
| sra-tools | 3.1.0 | Damona (Zenodo 10894351) |
| wget | 1.24.5 | Damona (Zenodo 14415567) |
| Trimmomatic | 0.39 | Damona (Zenodo 22145522) |
| Salmon | 1.2.0 | Damona (Zenodo 22144681) |
| gffread | 0.12.7 | Damona (Zenodo 8037039) |
