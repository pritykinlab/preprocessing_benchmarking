from utils import pipeline
from utils import processing_steps

import os
import scanpy as sc
import pandas as pd
import numpy as np
import scipy
from scipy.sparse import csr_matrix
from itertools import product
from sklearn.neighbors import NearestNeighbors




# Example pipeline configuration
test_params= [
    # Specify each step with its parameters
    (processing_steps.hvg, {'hvg_method': ['seurat', 'seurat_v3', 'pearson_residuals', 'min_cells', 'total_counts'], 'num_hvg': [2000, 3000, 4000, 5000, 6000, 7000, 8000, 12000, 16000]}),
    (processing_steps.norm, {'norm_method': ['log_zscore', 'pearson_residuals']}),
    (processing_steps.pca, {'max_pcs': [200]}),
    (processing_steps.evaluate, {'label_col': ['perturbation'],'num_nn': [20],'num_pcs_list': [10, 25, 30, 40, 50, 75, 100, 125, 150, 200]})
]

default_slurm_params = {
    'cpus-per-task': 4,
    'mem-per-cpu': '32GB',
}


# Uncomment the following line to run the pipeline with your specific input file and desired output directory
pipeline.run_pipeline(input_adata_file="/Genomics/pritykinlab/dillon/perturbseq/data/10_gene_analysis_v2/benchmark/for_pipeline.h5ad",
                            output_dir="/Genomics/pritykinlab/dillon/perturbseq/data/10_gene_analysis_v2/benchmark/results",
                            default_slurm_params = default_slurm_params,
                            pipeline_params=test_params,
                            verbose=True)
