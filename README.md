# scRNA-seq Analysis of Circulating Myeloid Cells in Neonatal HIE

**Master's Thesis (TFM) — Máster Universitario en Análisis de Datos Ómicos y Biología de Sistemas**  
Universidad de Sevilla · Universidad Internacional de Andalucía  
**Author:** Celia Galiano Manzanares

---

## Project Overview

This repository contains the full bioinformatics pipeline for a single-cell RNA sequencing (scRNA-seq) analysis of circulating myeloid cell populations in neonates with hypoxic-ischemic encephalopathy (HIE) compared to healthy controls.

The central biological question is whether myeloid subpopulations — classical monocytes, non-classical monocytes, intermediate monocytes, myeloid dendritic cells (mDCs), and plasmacytoid dendritic cells (pDCs) — undergo transcriptional reprogramming in HIE at 24–48 hours post-birth, a timepoint associated with the anti-inflammatory/immunoparalysis phase of the immune response rather than acute inflammation.

---

## Dataset

| Sample | Group | Notes |
|--------|-------|-------|
| C1_B3 | Control | |
| C2_B8 | Control | |
| HI1_B1 | HIE | |
| HI3_B2 | HIE | |
| HI4_B4 | HIE | |
| HI8_B5 | HIE | |
| HI9_B6 | HIE | Biological outlier: fetal hemoglobin overexpression (HBG1/HBG2), likely very premature neonate |
| HI11_B7 | HIE | **Excluded post-QC** (insufficient cell recovery) |

**Final integrated dataset:** 41,135 cells across 30 clusters.  
**Myeloid subset:** ~9,673 cells across 10 subclusters.

All downstream analyses (DGE, GO enrichment) were run in parallel in two sensitivity versions:
- `con_H9` — including sample HI9_B6
- `sin_H9` — excluding sample HI9_B6 (**final submitted analyses use `sin_H9` only**)

---

## Repository Structure

```
.
├── 1-Preprocesamiento_scRNAseq_v4.Rmd          # Per-sample QC, SoupX, doublet removal
├── 2-Integracion_Harmony_scRNAseq.Rmd          # Multi-sample integration with Harmony
├── 3-Visualizacion_integracion_scRNAseq.Rmd    # Global UMAP, annotation, composition plots
├── 4_Integracion_Harmony_Monocitos.Rmd         # Myeloid subset re-clustering with Harmony
├── 5_Visualizacion_Caracterizacion_Monocitos_v2.Rmd  # Myeloid characterization, manual annotation
├── 6_DGE_monocitos_v3.Rmd                      # Differential gene expression (MAST) + GO enrichment
├── 7_Pseudobulk_Monocitos_sinH9.Rmd            # Exploratory pseudobulk (distance/PCA)
└── README.md
```

---

## Pipeline Overview

The analysis is structured as seven sequential scripts:

### Script 1 — Per-sample Preprocessing
**`1-Preprocesamiento_scRNAseq_v4.Rmd`**

Executed independently for each of the 8 samples (one sample per run, controlled by the `nombre_muestra` parameter in the `init` chunk).

Steps:
1. **Ambient RNA correction** with SoupX (filtered + raw CellRanger matrices required). Estimated contamination fraction (ρ) reported per sample.
2. **Seurat object creation** with `min.cells = 3`, `min.features = 200`.
3. **Quality control filtering:** `200 < nFeature_RNA < 5000` and `percent.MT < 5%`.
4. **Erythrocyte removal** based on normalized HBB expression (threshold: `HBB ≥ 1`).
5. **Standard Seurat workflow:** normalization (LogNormalize, scale factor 10⁴), HVG selection (VST, 2000 features), scaling, PCA (dims 1:15).
6. **Doublet detection and removal** with scDblFinder (`set.seed(123)`).
7. **Automated cell type annotation** with SingleR using Monaco Immune Data (`label.fine`) and Human Primary Cell Atlas (`label.main`) as references.
8. Output: one `*_final.rds` Seurat object per sample.

> HI11_B7 was excluded after this step due to insufficient cell recovery.

---

### Script 2 — Multi-sample Integration
**`2-Integracion_Harmony_scRNAseq.Rmd`**

Integrates the 7 per-sample objects into a single dataset.

Steps:
1. Merge of 7 Seurat objects (`merge()` with `add.cell.ids`).
2. Re-normalization and PCA on the merged object (2000 HVGs, dims 1:20).
3. **Harmony integration** (`RunHarmony()`, `group.by.vars = "orig.ident"`, theta = 2, lambda = 1, max_iter = 10, dims = 1:20).
4. UMAP, neighbor graph, and clustering (resolution 0.8) on Harmony embeddings.
5. Second Harmony run correcting for both sample and condition (diagnostic only, stored as `harmony_doble`).
6. SingleR re-annotation on the integrated object (Monaco `label.fine` and `label.main`).
7. Myeloid subset extraction with double filter: clusters 1, 3, 9, 16, 24, 25, 26, 27 AND Monaco Main label in `{Monocytes, Macrophages, Dendritic cells, Microglia, Progenitors}`.
8. Output: `seu_integrado_harmony.rds`, `seu_mieloide.rds`.

---

### Script 3 — Global Visualization
**`3-Visualizacion_integracion_scRNAseq.Rmd`**

Generates publication-quality figures for the global integrated dataset.

Outputs: QC violin plots, UMAP panels (clusters, sample, condition), cell type composition barplots (per sample and Control vs HIE), DEG volcano plot (global Wilcoxon), heatmaps (markers per cluster, DEGs per cell type, cluster-vs-type composition).

---

### Script 4 — Myeloid Subset Re-clustering
**`4_Integracion_Harmony_Monocitos.Rmd`**

Dedicated Harmony re-integration of the myeloid subset to resolve monocyte subpopulations.

Steps:
1. Re-normalization, HVG selection (2000), PCA (dims 1:15) on the myeloid subset.
2. Harmony integration by sample (theta = 2, lambda = 1, dims 1:15).
3. Clustree analysis across resolutions 0.1–0.8 → final resolution: **0.4** (10 clusters).
4. Canonical marker visualization (CD14, FCGR3A, CX3CR1, S100A8/9, HLA-DRA, SELL, FCGR1A, CCR2) via FeaturePlots, VlnPlots, DotPlots.
5. SingleR re-annotation on the myeloid subset.
6. Output: `seu_monocitos_integrado.rds`.

---

### Script 5 — Myeloid Characterization and Manual Annotation
**`5_Visualizacion_Caracterizacion_Monocitos_v2.Rmd`**

Manual cluster identity assignment based on canonical marker expression.

| Cluster(s) | Assigned Identity |
|------------|-------------------|
| 0, 1, 4, 8 | Classical monocytes |
| 2 | Non-classical monocytes |
| 3, 5 | Intermediate monocytes |
| 6 | Myeloid DCs |
| 9 | Plasmacytoid DCs |
| 7 | Platelet mix |

> **Note on cluster 2 (Non-classical monocytes):** Manually reassigned despite 94.2% automated SingleR label as monocytes. Assignment based on marker inspection (FCGR3A, CX3CR1, CDKN1C) and topological UMAP isolation.

Outputs: annotated UMAP, composition barplots per sample and condition, DotPlots of identity markers (with/without condition split), Excel table of cluster identities.

---

### Script 6 — Differential Gene Expression and GO Enrichment
**`6_DGE_monocitos_v3.Rmd`**

DGE and functional enrichment analysis, run in parallel for `con_H9` and `sin_H9` scenarios.

**DGE parameters:**
- Test: **MAST** (`FindMarkers()`)
- Comparison: HIE vs Control within each myeloid subpopulation
- `min.pct = 0.01`, `logfc.threshold = 1` (log2FC)
- Multiple testing correction: BH (FDR)
- Minimum group size: 10 cells per condition

> Intermediate monocytes (n = 2 Control cells) and plasmacytoid DCs (0 Control cells) did not meet the minimum threshold and were excluded from DGE.

**GO enrichment parameters:**
- `clusterProfiler::enrichGO()`, ontology = BP, BH correction, FDR < 0.05
- Background: all genes tested per cell type (universe)
- Semantic reduction: `rrvgo` (threshold = 0.7, method = "Rel")
- Outputs: dotplots, treemaps (PDF), Excel with full and reduced GO terms

**Final submitted results:** `sin_H9` scenario only.

---

### Script 7 — Exploratory Pseudobulk
**`7_Pseudobulk_Monocitos_sinH9.Rmd`**

Exploratory pseudobulk analysis of the myeloid compartment as a whole (no subtype stratification).

Steps:
1. Count aggregation per sample: `AggregateExpression()` on `orig.ident`.
2. DESeq2 object construction and VST normalization (`blind = TRUE`).
3. Euclidean distance heatmap between samples.
4. Pearson correlation heatmap between samples.
5. PCA on VST-transformed pseudobulk profiles.

> **Important:** This analysis is **exploratory only**. DESeq2 inferential testing (Wald test) was not applied because the study design has n = 2 controls, which is below the minimum required for valid pseudobulk differential expression (≥ 3 replicates per group). The pseudobulk framework is correctly classified as descriptive (sample-level QC and condition separation assessment).

---

## Software and Versions

| Package | Version |
|---------|---------|
| R | 4.4.2 |
| Seurat | 5.4.0 |
| SeuratObject | 5.3.0 |
| SoupX | 1.6.2 |
| scDblFinder | 1.20.2 |
| SingleR | 2.8.0 |
| harmony | 1.2.0 |
| MAST | — |
| clusterProfiler | 4.x |
| rrvgo | — |
| DESeq2 | — |
| clustree | 0.5.1 |
| openxlsx | — |
| ggplot2 | — |
| patchwork | — |

**Reference datasets:**
- Monaco Immune Data (Monaco et al., *Cell Reports*, 2019; GSE107011) — primary annotation
- Human Primary Cell Atlas — validation annotation

---

## Cell Type Color Palette

| Cell Type | Hex Color |
|-----------|-----------|
| Classical monocytes | `#CAFF70` |
| Non-classical monocytes | `#97FFFF` |
| Intermediate monocytes | `#FFEC8B` |
| Myeloid DCs | `#EE7AE9` |
| Plasmacytoid DCs | `#CDBA96` |
| Platelet mix | `#FFAEB9` |

---

## Key Seurat Objects

| File | Contents |
|------|----------|
| `*_final.rds` (×7) | Per-sample preprocessed objects |
| `seu_integrado_harmony.rds` | Full 7-sample integrated object (41,135 cells, 30 clusters) |
| `seu_monocitos_integrado.rds` | Myeloid subset, annotated (con_H9 version) |
| `seu_mono_final.rds` | Myeloid subset, sin_H9 version (used for DGE/GO) |

---

## Input Data

Raw CellRanger output (filtered and raw feature-barcode matrices) required per sample:
```
samples/
└── {SAMPLE_NAME}/
    ├── filtered_feature_bc_matrix/
    └── raw_feature_bc_matrix/
```

Data can be accessed on request to the authors.

---

## Statistical Notes

- **n = 2 controls** is a known limitation of this dataset. All DGE results should be interpreted with caution and treated as exploratory.
- **No pseudobulk inferential framework** (DESeq2 Wald test) was applied due to the insufficient number of replicates per group.
- **MAST** at the single-cell level was used as the primary DGE method (Script 6).
- The `sin_H9` sensitivity analysis was conducted to assess the influence of a biological outlier (HI9_B6: excessive fetal hemoglobin expression, likely very premature neonate).

---

## References

- Stuart et al. (2019). Comprehensive integration of single-cell data. *Cell*. [Seurat]
- Korsunsky et al. (2019). Fast, sensitive and accurate integration of single-cell data with Harmony. *Nature Methods*.
- Aran et al. (2019). Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage. *Nature Immunology*. [SingleR]
- Monaco et al. (2019). RNA-seq signatures normalized by mRNA abundance allow absolute deconvolution of human cell types. *Cell Reports*.
- Young & Behjati (2020). SoupX removes ambient RNA contamination from droplet-based single-cell RNA sequencing data. *GigaScience*. [SoupX]
- Germain et al. (2021). Doublet identification in single-cell sequencing data using scDblFinder. *F1000Research*.
- Finak et al. (2015). MAST: a flexible statistical framework for assessing transcriptional changes and characterizing heterogeneity in single-cell RNA sequencing data. *Genome Biology*.
- Wu et al. (2021). clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. *Innovation*.
- Sayols (2020). rrvgo: a Bioconductor package for interpreting lists of Gene Ontology terms. *microPublication Biology*.
- Squair et al. (2021). Confronting false discoveries in single-cell differential expression. *Nature Communications*.
- Chakkarapani et al. (2025). [HIE clinical context].
- Reyes et al. (2020). An immune-cell signature of bacterial sepsis. *Nature Medicine*.

---

## License

This code is provided for academic reproducibility purposes. Please cite the corresponding thesis if you use or adapt any part of this pipeline.

