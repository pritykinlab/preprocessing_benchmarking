# Dillon To Do:
- normalize the metric
- many dataset from large dataset

# Yujie To Do:
- HVGs = [250, 500, 1000, 2000, 4000, 8000, max_genes]
- aggregate the results by doing zscore then plotting all of the results on one dataset (num_PCs x num_k x num_norm plots)

# TO DO 1
- identify the control cells
- where are the replicates are there any
- remove controls from the metric

# TO DO
- add sctransform
- highly variable genes need to go to the max value (edit the metadata to say how many genes are at max then pass that as the argument)

# Long Term To Do:
- subsets of the data (in rare subpopulations, I suspect that more HVG better because, a gene might not be HVG unless if you have enough cells of the small subpopulation)
