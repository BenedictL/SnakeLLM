# Workflow Pattern: Single-Cell RNA-seq Analysis
Category: scrna-seq
Analysis type: single cell, cell clustering, cell type annotation
Keywords: scRNA-seq, single cell, 10x Genomics, Cell Ranger, Seurat, Scanpy, clustering, UMAP, cell types

## Purpose
Profile gene expression at single-cell resolution. Identify distinct cell populations,
marker genes, and cell-type-specific expression patterns.

## Pipeline Steps (in order)
1. cellranger_count    — Cell Ranger: alignment + UMI counting (10x Genomics data)
   OR star_solo        — STARsolo: open-source alternative to Cell Ranger
2. quality_control     — Seurat or Scanpy: filter low-quality cells (nCount, nFeature, % mito)
3. normalization       — Normalize + log-transform count matrix
4. feature_selection   — Select highly variable genes (HVGs)
5. dimensionality_red  — PCA → UMAP/tSNE for visualization
6. clustering          — Leiden or Louvain graph clustering
7. marker_genes        — FindMarkers / rank_genes_groups: identify cluster markers
8. cell_annotation     — Annotate clusters using known markers or automated tools (SingleR, CellTypist)
9. diff_expression     — Pseudobulk DE between conditions (DESeq2 on aggregated counts)

## Key Decision Points
- Cell Ranger vs STARsolo: Cell Ranger is the 10x gold standard; STARsolo is free and nearly equivalent
- Seurat (R) vs Scanpy (Python): Seurat more established; Scanpy integrates better with Python ML ecosystem
- Clustering resolution: test multiple resolutions (0.2 to 1.5); validate with known markers