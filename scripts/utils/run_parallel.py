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

def run(row):
    output_base_dir = "/Genomics/pritykinlab/yujie/preprocessing_benchmarking/results"
    input_datasets_dir = "/Genomics/pritykinlab/yujie/preprocessing_benchmarking/datasets/harmonized_perturb_datasets"
    size = row['Size (GB)']
    dataset_file = row['Dataset']
    # Construct the full path to the dataset file
    dataset_path = os.path.join(input_datasets_dir, dataset_file)

    # Generate a unique output directory for each dataset
    dataset_name = dataset_file[:-5]  # Remove '.h5ad' extension for folder name
    output_dir = os.path.join(output_base_dir, dataset_name)

    # Create the output directory if it does not exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    test_params= [
        # Specify each step with its parameters
        (processing_steps.hvg, {'hvg_method': ['min_cells'], 'num_hvg': [1000, 8000]}),
        (processing_steps.norm, {'norm_method': ['log_zscore', 'pearson_residuals']}),
        (processing_steps.pca, {'max_pcs': [500]}),
        (processing_steps.evaluate, {'label_col': ['perturbation'],'num_nn': [20],'num_pcs_list': [10, 500]})
    ]

    default_slurm_params = {
        'cpus-per-task': 2,
        'mem-per-cpu': '32G',
    }
    print(default_slurm_params)

    # Run the pipeline with the current dataset file and unique output directory
    pipeline.run_pipeline(input_adata_file=dataset_path,
                          output_dir=output_dir,
                          default_slurm_params=default_slurm_params,
                          pipeline_params=test_params,
                          verbose=True,
                          parallel_type='regular')
