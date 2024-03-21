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

def run(dataset, default_slurm_params):
    name = dataset.split("/")[-1].split(".")[0]
    output_base_dir = f"/Genomics/pritykinlab/yujie/preprocessing_benchmarking/results/harmonized_perturb_seq_new_pipeline2/{name}"

    test_params= [
        # Specify each step with its parameters
        (processing_steps.clean, {}),
        (processing_steps.hvg_norm, {'hvg_norm_combo': ['Pearson Residual + Pearson Residual', 'Pearson Residual + norm_log_zscore', 'seurat + norm_log_zscore', 'seurat + norm_log'], 'num_hvg': [150, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000]}),
        (processing_steps.pca, {'max_pcs': [100]}),
        (processing_steps.evaluate, {'label_col': ['perturbation'],'num_nn': [20, 40],'num_pcs': [25, 50, 100]})
    ]

    print(default_slurm_params)

    # Run the pipeline with the current dataset file and unique output directory
    pipeline.run_pipeline(input_adata_file=dataset,
                          output_dir=output_base_dir,
                          default_slurm_params=default_slurm_params,
                          pipeline_params=test_params,
                          verbose=True,
                          remove_intermediate=True,
                          parallel_type='slurm')


def run_10x_CITE_seq(dataset, default_slurm_params):
    name = dataset.split("/")[-1].split(".")[0]
    output_base_dir = f"/Genomics/pritykinlab/yujie/preprocessing_benchmarking/results/10x_CITE_seq/{name}"
    protein_h5ad_file = f"/Genomics/pritykinlab/yujie/preprocessing_benchmarking/results/10x_CITE_seq/{name}/intermediate_files/protein.h5ad"

    test_params= [
        # Specify each step with its parameters
        (processing_steps.clean_CITE_seq, {'protein_processed_out_file': protein_h5ad_file}),
        (processing_steps.hvg_norm, {'hvg_norm_combo': ['Pearson Residual + Pearson Residual', 'Pearson Residual + norm_log_zscore', 'seurat + norm_log_zscore', 'seurat + norm_log'], 'num_hvg': [150, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000]}),
        (processing_steps.pca, {'max_pcs': [100]}),
        (processing_steps.evaluate_CITE_seq, {'num_nn': [20, 40],'num_pcs': [25, 50, 100], 'protein_h5ad_file': protein_h5ad_file})
    ]

    print(default_slurm_params)

    # Run the pipeline with the current dataset file and unique output directory
    pipeline.run_pipeline(input_adata_file=dataset,
                          output_dir=output_base_dir,
                          default_slurm_params=default_slurm_params,
                          pipeline_params=test_params,
                          verbose=True,
                          remove_intermediate=True,
                          parallel_type='slurm')


