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

output_base_dir = "/Genomics/pritykinlab/yujie/preprocessing_benchmarking/results"
input_datasets_dir = "/Genomics/pritykinlab/yujie/preprocessing_benchmarking/datasets/harmonized_perturb_datasets"

metadata_df = pd.read_csv("../dataset_metadata/selected_datasets_summary.tsv", sep='\t')
metadata_df = metadata_df.sort_values(['Size (GB)'])

# Loop through each file in the input datasets directory
for idx, row in metadata_df.iterrows():
    size = row['Size (GB)']
    dataset_file = row['Dataset']
    # Construct the full path to the dataset file
    dataset_path = os.path.join(input_datasets_dir, dataset_file)
    
    # Skip if it's not a file or if it doesn't end with '.h5ad'
    if not os.path.isfile(dataset_path) or not dataset_file.endswith(".h5ad"):
        continue
    
    # Generate a unique output directory for each dataset
    dataset_name = dataset_file[:-5]  # Remove '.h5ad' extension for folder name
    output_dir = os.path.join(output_base_dir, dataset_name)
    
    # Create the output directory if it does not exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Example pipeline configuration
#    test_params= [
#        # Specify each step with its parameters
#        (processing_steps.hvg, {'hvg_method': ['seurat', 'cell_ranger', 'pearson_residuals', 'min_cells', 'total_counts'], 'num_hvg': [1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000]}),
#        (processing_steps.norm, {'norm_method': ['log_zscore', 'pearson_residuals']}),
#        (processing_steps.pca, {'max_pcs': [500]}),
#        (processing_steps.evaluate, {'label_col': ['perturbation'],'num_nn': [20],'num_pcs_list': [10, 25, 50, 75, 100, 125, 150, 175, 200, 350, 400, 500]})
#    ]
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



