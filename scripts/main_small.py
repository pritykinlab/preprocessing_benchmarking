import os
import scanpy as sc
import pandas as pd
import numpy as np
import scipy
from scipy.sparse import csr_matrix
from itertools import product
from sklearn.neighbors import NearestNeighbors
from utils import pipeline
from utils import processing_steps

os.environ["OPENBLAS_NUM_THREADS"] = "64"

output_base_dir = "/Genomics/pritykinlab/yujie/preprocessing_benchmarking/test"
input_datasets_dir = "/Genomics/pritykinlab/yujie/preprocessing_benchmarking/datasets/harmonized_perturb_datasets"

# Function to get the size of the dataset
def get_dataset_size(file_path):
    return os.path.getsize(file_path)

# Find the smallest dataset
smallest_size = None
smallest_dataset_file = None

for file in os.listdir(input_datasets_dir):
    if file.endswith(".h5ad"):
        file_path = os.path.join(input_datasets_dir, file)
        dataset_size = get_dataset_size(file_path)
        
        if smallest_size is None or dataset_size < smallest_size:
            smallest_size = dataset_size
            smallest_dataset_file = file

# Continue with the pipeline for the smallest dataset
if smallest_dataset_file:
    dataset_path = os.path.join(input_datasets_dir, smallest_dataset_file)
    dataset_name = smallest_dataset_file[:-5]  # Remove '.h5ad' extension for folder name
    output_dir = os.path.join(output_base_dir, dataset_name)
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    test_params= [
        # Specify each step with its parameters
        (processing_steps.hvg_norm, {'hvg_norm_combo': ['Pearson Residual + Pearson Residual', 'Pearson Residual + log_zscore', 'seurat + log_zscore', 'seurat + log'], 'num_hvg': [1000, 5000]}),
        (processing_steps.pca, {'max_pcs': [100]}),
        (processing_steps.evaluate, {'label_col': ['perturbation'],'num_nn': [20, 40],'num_pcs_list': [25, 100]})
    ]

    default_slurm_params = {
        'cpus-per-task': 2,
        'mem-per-cpu': '32G',
    }
    
    # Run the pipeline with the specified dataset file and unique output directory
    pipeline.run_pipeline(input_adata_file=dataset_path,
                          output_dir=output_dir,
                          default_slurm_params=default_slurm_params,
                          pipeline_params=test_params,
                          verbose=True,
                          remove_intermediate=False,
                          parallel_type='regular')
else:
    print("No .h5ad files found in the directory.")
