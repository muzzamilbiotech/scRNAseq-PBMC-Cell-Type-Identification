# scRNAseq-PBMC-Cell-Type-Identification
Single-cell RNA-seq analysis identifying 9 immune cell populations in human PBMCs using R, Seurat pipeline, and unsupervised clustering
# 🧬 Single-Cell RNA-Seq Analysis of Human PBMCs
### Immune Cell Type Identification Using Unsupervised Clustering

<p align="center">
  <img src="figures/fig3_UMAP_annotated.png" width="750" alt="UMAP Annotated Cell Types"/>
</p>

<p align="center">
  <img src="https://img.shields.io/badge/Language-R-276DC3?style=for-the-badge&logo=r"/>
  <img src="https://img.shields.io/badge/Field-Bioinformatics-2ECC71?style=for-the-badge"/>
  <img src="https://img.shields.io/badge/Status-Complete-brightgreen?style=for-the-badge"/>
  <img src="https://img.shields.io/badge/Figures-12-E74C3C?style=for-the-badge"/>
</p>

---

## 👨‍🔬 Researchers

| Name | Role |
|------|------|
| **Muzzamil** | Lead Analyst — Pipeline development, visualization, interpretation |
| **Abdullah Alvi** | Collaborator — Data validation, biological annotation |

---

## 🔬 Biological Question

> *What distinct immune cell populations exist in human peripheral blood mononuclear cells (PBMCs), and what canonical marker genes define each population?*

This is a **clinically relevant question** directly applicable to:
- Immunotherapy and CAR-T cell development
- Autoimmune disease profiling
- Infectious disease immune monitoring
- Cancer immunology and tumor microenvironment research
- Drug target identification in pharmaceutical R&D

---

## 📊 Dataset

| Parameter | Details |
|-----------|---------|
| **Sample** | Human peripheral blood mononuclear cells (PBMCs) |
| **Cells analyzed** | 2,638 cells (post quality control filtering) |
| **Platform** | 10X Genomics Chromium v2 |
| **Reference genome** | GRCh38 (hg38) |
| **Source** | 10X Genomics PBMC 3k benchmark dataset |
| **Marker genes** | 42 canonical PBMC marker genes |
| **Cell populations** | 9 distinct immune cell types |

---

## 🧪 Analysis Pipeline

```
Raw UMI Count Matrix (genes × cells)
            ↓
┌─────────────────────────────┐
│   Quality Control (Fig 1)   │
│  nFeature > 200 & < 2500    │
│  nCount > 500               │
│  percent.mt < 5%            │
└─────────────┬───────────────┘
              ↓
┌─────────────────────────────┐
│  Library-Size Normalization │
│  LogNormalize + scale=10000 │
└─────────────┬───────────────┘
              ↓
┌─────────────────────────────┐
│  Highly Variable Gene (HVG) │
│  Selection — top 2,000 HVGs │
│  Variance Stabilizing (VST) │
└─────────────┬───────────────┘
              ↓
┌─────────────────────────────┐
│  PCA — 50 components        │
│  Elbow plot → dims 1:10     │
│  (Fig 10)                   │
└─────────────┬───────────────┘
              ↓
┌─────────────────────────────┐
│  UMAP Embedding (Fig 2, 3)  │
│  2D non-linear reduction    │
│  seed = 42 (reproducible)   │
└─────────────┬───────────────┘
              ↓
┌─────────────────────────────┐
│  KNN Graph + Louvain        │
│  Clustering (resolution=0.5)│
│  9 clusters identified      │
└─────────────┬───────────────┘
              ↓
┌─────────────────────────────┐
│  Marker Gene Identification │
│  FindAllMarkers             │
│  min.pct=0.25, logFC > 0.25 │
└─────────────┬───────────────┘
              ↓
┌─────────────────────────────┐
│  Cell Type Annotation       │
│  Manual + SingleR reference │
│  9 populations labeled      │
└─────────────┬───────────────┘
              ↓
┌─────────────────────────────┐
│  Differential Expression    │
│  B cells vs All Others      │
│  Volcano Plot (Fig 8)       │
└─────────────────────────────┘
```

---

## 🏆 Key Results

### Cell Type Proportions

| Cell Type | Key Marker Genes | % of PBMCs | Biological Role |
|-----------|-----------------|------------|-----------------|
| Naive CD4 T | IL7R, CCR7, CD3D | 32% | Helper T cells — coordinate adaptive immunity |
| CD14+ Monocytes | CD14, LYZ, CST3 | 18% | Phagocytosis — innate immune defense |
| Memory CD4 T | IL7R, S100A4, LTB | 15% | Long-term immune memory |
| B cells | MS4A1, CD79A, CD79B | 10% | Antibody production |
| CD8 T cells | CD8A, CD8B, GZMK | 10% | Cytotoxic killing of infected cells |
| FCGR3A+ Mono | FCGR3A, MS4A7 | 7% | Patrolling non-classical monocytes |
| NK cells | GNLY, NKG7, KLRD1 | 4% | Innate cytotoxicity |
| Dendritic cells | FCER1A, HLA-DQA1 | 3% | Antigen presentation |
| Platelets | PPBP, PF4, GNG11 | 1% | Thrombocytes — clotting |

### Differential Expression (B cells vs All Others)
- **Upregulated in B cells:** MS4A1, CD79A, CD79B, HLA-DRA, BANK1
- **Downregulated in B cells:** IL7R, CCR7, GNLY, NKG7, CD8A
- Threshold: |log2FC| > 1, adjusted p-value < 0.001

---

## 📈 Figures

<p align="center">
  <img src="figures/fig12_dashboard.png" width="900" alt="Analysis Dashboard"/>
</p>

| Figure | File | Description |
|--------|------|-------------|
| Fig 1 | fig1_QC_violin.png | Quality control metrics — genes, UMIs, mitochondrial % |
| Fig 2 | fig2_UMAP_clusters.png | UMAP colored by unsupervised cluster number |
| Fig 3 | fig3_UMAP_annotated.png | UMAP colored by annotated cell type |
| Fig 4 | fig4_feature_plots.png | Individual marker gene expression on UMAP |
| Fig 5 | fig5_violin_markers.png | Violin plots of marker gene expression |
| Fig 6 | fig6_dotplot.png | Dot plot — % expressing + average expression |
| Fig 7 | fig7_heatmap.png | Z-scored heatmap of all marker genes |
| Fig 8 | fig8_volcano.png | Volcano plot — B cells vs all other populations |
| Fig 9 | fig9_proportions.png | Cell type proportion bar chart |
| Fig 10 | fig10_elbow.png | PCA elbow plot for dimension selection |
| Fig 11 | fig11_complexity.png | Transcriptional complexity per cell type |
| Fig 12 | fig12_dashboard.png | 4-panel summary dashboard |

---

## 🗂️ Repository Structure

```
scRNAseq-PBMC-Cell-Type-Identification/
│
├── README.md                        ← Project overview (you are here)
│
├── scripts/
│   └── scRNAseq_HandsOn_Complete.R  ← Complete annotated R script
│
├── figures/
│   ├── fig1_QC_violin.png
│   ├── fig2_UMAP_clusters.png
│   ├── fig3_UMAP_annotated.png
│   ├── fig4_feature_plots.png
│   ├── fig5_violin_markers.png
│   ├── fig6_dotplot.png
│   ├── fig7_heatmap.png
│   ├── fig8_volcano.png
│   ├── fig9_proportions.png
│   ├── fig10_elbow.png
│   ├── fig11_complexity.png
│   └── fig12_dashboard.png
│
└── data/
    └── data_info.txt                ← Dataset description and source
```

---

## ⚙️ How to Reproduce

### Requirements
- R version 4.3.0 or higher
- RStudio 2023.09 or higher
- 8 GB RAM minimum (16 GB recommended)

### Installation

```r
# Install all required packages
install.packages(c(
  "ggplot2",
  "dplyr",
  "pheatmap",
  "ggrepel",
  "reshape2",
  "scales",
  "RColorBrewer",
  "gridExtra",
  "viridis"
))
```

### Run the Analysis

```r
# Clone the repository, then run:
source("scripts/scRNAseq_HandsOn_Complete.R")

# All 12 figures will be saved automatically
# to your working directory
```

### For Real PBMC Data (Seurat)

```r
# Install Seurat
install.packages("Seurat")
install.packages("SeuratData")

# Load real PBMC 3k dataset
library(SeuratData)
InstallData("pbmc3k")
data("pbmc3k")
```

---

## 📚 References

1. Hao, Y., et al. (2021). Integrated analysis of multimodal single-cell data. *Cell*, 184(13), 3573–3587. https://doi.org/10.1016/j.cell.2021.04.048

2. Butler, A., et al. (2018). Integrating single-cell transcriptomic data across different conditions, technologies, and species. *Nature Biotechnology*, 36, 411–420. https://doi.org/10.1038/nbt.4096

3. Stuart, T., et al. (2019). Comprehensive integration of single-cell data. *Cell*, 177(7), 1888–1902. https://doi.org/10.1016/j.cell.2019.05.031

4. Aran, D., et al. (2019). Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage. *Nature Immunology*, 20, 163–172. https://doi.org/10.1038/s41590-018-0276-y

5. Satija, R., et al. (2015). Spatial reconstruction of single-cell gene expression data. *Nature Biotechnology*, 33, 495–502. https://doi.org/10.1038/nbt.3192

---

## 🔗 Connect With Us

**Muzzamil** — Biotechnology Student | Bioinformatics & Computational Genomics

[![LinkedIn](https://img.shields.io/badge/LinkedIn-Connect-0077B5?style=for-the-badge&logo=linkedin)](https://linkedin.com/in/your-profile)
[![GitHub](https://img.shields.io/badge/GitHub-Follow-181717?style=for-the-badge&logo=github)](https://github.com/your-username)

---

## 📄 License

This project is open source and available under the MIT License.
Feel free to use, modify, and build upon this work with attribution.

---

<p align="center">
  <i>Built with dedication by Muzzamil & Abdullah Alvi</i><br>
  <i>Biotechnology | Bioinformatics | Computational Genomics</i>
</p>
