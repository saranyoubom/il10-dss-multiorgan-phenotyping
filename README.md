# Multi-organ metabolic and transcriptional phenotyping of IL-10-deficient mice under graded intestinal inflammatory challenge

Data and analysis code supporting the manuscript submitted to *Scientific Reports*
(Submission ID: 62ee1d26-f250-4408-aa6a-409aeb64f63d).

## Study groups

- **CTRL** — C57BL/6J background control mice
- **IL10-KO** — untreated IL-10-knockout mice (B6.129P2-Il10tm1Cgn)
- **DSS-IL10-KO** — IL-10-knockout mice receiving recurrent DSS treatment (2 cycles, weeks 0 and 6, with a 5-week remitting interval)

n = 5 animals per group for whole-animal physiological endpoints; n = 4 per group for RT-qPCR
gene expression (after exclusion of amplification-QC failures). Male mice only (see manuscript
Methods for rationale).

## Folder structure

| Folder | Content |
|---|---|
| `Fig1_Intestinal/` | Body weight, FITC-dextran intestinal permeability, albumin/globulin ratio, large-intestine (colon) RT-qPCR gene expression, colon gene-pair correlation matrix |
| `Fig2_Pancreas/` | Fasting blood glucose, fasting C-peptide, FCP/FBG ratio, pancreatic beta-cell/glucose-handling/T2D-susceptibility RT-qPCR gene expression, pancreas gene-pair correlation matrix |
| `Fig3_Immune/` | 40-plex serum cytokine array (raw + log2, hierarchical clustering), serum blood chemistry (validated against Otto et al. 2016 reference ranges), splenic RT-qPCR gene expression, spleen gene-pair correlation matrix |
| `Fig4_Organs/` | Hepatic and renal RT-qPCR gene expression, liver and kidney gene-pair correlation matrices |
| `Fig5_InterOrgan/` | Inter-organ (cross-tissue) gene-pair correlation analysis — see note below |
| `Supplementary/` | Blood chemistry validation detail (Supplementary Fig. 1) and an independent pancreas gene-panel correlation check (Supplementary Fig. 2) |

## Fig5_InterOrgan — two related analyses

- `MetabolicImmune_Axis_Correlation_Results.csv` / `_Significant_Correlations.csv` /
  `MetabolicImmune_Axis_Correlation_Analysis.r` — the full pairwise correlation sweep across
  all RT-qPCR genes measured in each organ (1,210 gene-pair tests across the 5 organs), pooled
  across all three groups. This is the source analysis behind the published Fig. 5 panels.
- `Rebuttal_Pooled_InterOrgan_Correlations_BHcorrected.csv` /
  `Rebuttal_Stratified_InterOrgan_Correlations_PerGroup.csv` — a curated reanalysis restricted
  to the specific gene panel discussed in the main text per organ (287 gene-pair tests),
  adding Benjamini–Hochberg correction across all tests and a within-group-stratified
  correlation pass. Produced during peer review in response to reviewer comments on multiple
  testing and pooled-vs-within-group correlation structure; reported in Supplementary Table S3.

## File naming conventions

| Suffix | Meaning |
|---|---|
| `_Data_*` | Raw/individual-animal measurements |
| `_Summary_*` | Group-level descriptive statistics (mean, SD, median, IQR, min/max) |
| `_Normality_*` / `_Normality_Tests` / `_Normality_Results` | Shapiro–Wilk normality test results, used to select parametric vs. non-parametric tests |
| `_Statistical_Tests` / `_Test_Decision(s)` | Primary hypothesis test results and the parametric/non-parametric test selection made per the normality result |
| `_Posthoc_*` | Post-hoc pairwise comparisons, Benjamini–Hochberg (BH) corrected where noted |
| `_FoldChange_*` | Fold-change values relative to CTRL |
| `_Correlation_Matrix_*` / `_PValue_Matrix_*` / `_N_Matrix_*` | Pearson correlation coefficient, p-value, and sample-size (N) matrices from `Hmisc::rcorr()`, exported separately per the analysis convention used throughout this study |
| `_Pairwise_Correlations_*` | The same correlation results in long (one row per gene pair) format |
| `.r` | R analysis script that generated the accompanying CSV(s) and the corresponding manuscript figure panel |

## Gene abbreviations

Ins1/Ins2 (Insulin 1/2), Mafa (MAF BZIP Transcription Factor), Glp1r (glucagon-like peptide 1
receptor), Ghrhr (growth hormone releasing hormone receptor), Insr (insulin receptor), Tbc1d4
(TBC1 domain family member 4), Tcf19 (transcription factor 19), Srr (serine racemase), Zfand6
(zinc finger AN1-type domain 6), Glut1/Slc2a1 (solute carrier family 2 member 1), Gapdh
(glyceraldehyde-3-phosphate dehydrogenase, reference gene).

## Statistical methods

Shapiro–Wilk normality testing determined parametric (paired/Welch t-test) vs. non-parametric
(Kruskal–Wallis/Wilcoxon) test selection per comparison. Post-hoc pairwise comparisons used
Benjamini–Hochberg correction. Pearson correlation matrices report r, p-value, and N separately.
Full statistical methods are described in the manuscript Methods section.

## Related manuscript

Data Availability section of the manuscript points to this deposit's DOI. Raw histology images
are available from the corresponding author on reasonable request (not included here).

## License

Released under CC BY 4.0 unless the corresponding author specifies otherwise at deposit time.

## Contact

Corresponding author: see manuscript author list.
