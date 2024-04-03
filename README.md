# Dillon To Do:
- normalize the metric
- remove controls from the metric

# Yujie To Do:
- run a test with max hvg, test num_PCs: make plots: 20, 30, 40, 50, 60, 70, 80. 90, 100, 125, 150, 175, 200
- many dataset from large dataset
- identify the control cells (add a new boolean column into adata.obs "is_control_or_not")
- prepare a slide to explain zscore methods
- automate the plot script
- add num_PCs [5, 500] to show baseline

# TO DO 1
- where are the replicates are there any


# Long Term To Do:
- subsets of the data (in rare subpopulations, I suspect that more HVG better because, a gene might not be HVG unless if you have enough cells of the small subpopulation)

# More Things to benchmark
## incorporate new lower dimensional embeddings
- https://www.biorxiv.org/content/biorxiv/early/2024/03/27/2024.03.23.586420.full.pdf (GLM PCA)
- L1 and L2 PCA
- SCVI

## Normalization Techniques
- Weinreb et al (2018) used a simple extension of CPM that excludes genes that account for at least 5% of the total counts in any cell, when calculating their size factors
- Scran

# Miscellaneous
- Ground truth is known pathways

# Atlas Papers
- https://www.nature.com/articles/s41588-024-01688-9
