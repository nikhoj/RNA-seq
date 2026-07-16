# RNA-seq analysis report - `airway`

_Generated 2026-07-14 | thresholds: padj < 0.05, |log2FC| >= 1.0_

## Differential expression highlights

- **dexamethasone_vs_untreated**: 936 significant (486 up, 450 down). Top genes: SPARCL1, SAMHD1, MAOA, GPX3, DUSP1, SERPINA3, CACNB2, VCAM1.

## Phase 1 - QC

### count_densities

![count_densities](results/figures/phase1_qc/count_densities.png)

### detected_genes

![detected_genes](results/figures/phase1_qc/detected_genes.png)

### library_size

![library_size](results/figures/phase1_qc/library_size.png)

## Phase 2 - EDA / normalization

### normalization_boxplots

![normalization_boxplots](results/figures/phase2_qc/normalization_boxplots.png)

### pca

![pca](results/figures/phase2_qc/pca.png)

### sample_distance_heatmap

![sample_distance_heatmap](results/figures/phase2_qc/sample_distance_heatmap.png)

## Phase 3 - Differential expression

### de_summary_bar

![de_summary_bar](results/figures/phase3_de/de_summary_bar.png)

### ma_dexamethasone_vs_untreated

![ma_dexamethasone_vs_untreated](results/figures/phase3_de/ma_dexamethasone_vs_untreated.png)

### top_genes_dexamethasone_vs_untreated

![top_genes_dexamethasone_vs_untreated](results/figures/phase3_de/top_genes_dexamethasone_vs_untreated.png)

### volcano_dexamethasone_vs_untreated

![volcano_dexamethasone_vs_untreated](results/figures/phase3_de/volcano_dexamethasone_vs_untreated.png)

## Phase 5 - Enrichment

### gsea_dexamethasone_vs_untreated_Adipogenesis

![gsea_dexamethasone_vs_untreated_Adipogenesis](results/figures/phase5_enrichment/gsea_dexamethasone_vs_untreated_Adipogenesis.png)

### gsea_dexamethasone_vs_untreated_Androgen_Response

![gsea_dexamethasone_vs_untreated_Androgen_Response](results/figures/phase5_enrichment/gsea_dexamethasone_vs_untreated_Androgen_Response.png)

### ora_dotplot_dexamethasone_vs_untreated

![ora_dotplot_dexamethasone_vs_untreated](results/figures/phase5_enrichment/ora_dotplot_dexamethasone_vs_untreated.png)

## Tables

- [airway_DE_dexamethasone_vs_untreated.csv](results/tables/airway_DE_dexamethasone_vs_untreated.csv)
- [airway_GSEA_dexamethasone_vs_untreated.csv](results/tables/airway_GSEA_dexamethasone_vs_untreated.csv)
- [airway_ORA_dexamethasone_vs_untreated.csv](results/tables/airway_ORA_dexamethasone_vs_untreated.csv)
- [airway_temporal_clusters.csv](results/tables/airway_temporal_clusters.csv)
- [airway_timecourse_LRT.csv](results/tables/airway_timecourse_LRT.csv)
