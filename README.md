# Benchmarking computational analysis of scRNA-seq data using external biological context
Repository for systematic benchmarking of computational methods and parameters for scRNA-seq data analysis.

## Directory:
- `analysis/` - Scripts to visualize trends and patterns from resulting metrics files. 
- `dataset_metadata/` - Metadata files containing information in datasets attributes, such as number of cells/genes, file size, and whether it has raw counts.
- `download_dataset_scripts/` - Notebooks for downloading datasets from original sources, as well as steps for preprocessing raw data into a format ready for benchmarking analysis.
  - `CITE_seq/` - Notebooks to download CITE-seq datasets from original publications and steps to generate processed adata files.
  - `TCR_BCR_seq/` - Notebooks to download scTCR-seq and scBCR-seq datasets from original publications and steps to generate processed adata files.
  - `lineage_tracing/` - Notebooks to download lineage tracing datasets from original publications and steps to generate processed adata files.
  - `perturbseq/` - Notebooks to download Perturb-seq datasets from original publications and steps to generate processed adata files.
- `pilot_data/` - Notebooks for the analysis and generation of Figures 1 and 3.
- `scripts/` - Scripts for various downstream analysis.
  - `misc/` - Miscellaneous scripts such as preliminary data analysis.
  - `utils/` - Including functions used in the benchmarking pipeline as well as scripts for generating results in parallel.

