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
from utils import run_parallel

output_base_dir = "/Genomics/pritykinlab/yujie/preprocessing_benchmarking/results/harmonized_perturb_seq_test"
input_datasets_dir = "/Genomics/pritykinlab/yujie/preprocessing_benchmarking/datasets/harmonized_perturb_datasets"

metadata_df = pd.read_csv("../dataset_metadata/selected_datasets_summary.tsv", sep='\t')

# Exclude the dataset "TianKampmann2019_iPSC" from the DataFrame
metadata_df = metadata_df[metadata_df['Dataset'] != 'TianKampmann2019_iPSC.h5ad']

# metadata_df = metadata_df[metadata_df['Dataset'] != 'ReplogleWeissman2022_rpe1.h5ad']

# metadata_df = metadata_df[metadata_df['Dataset'] != 'FrangiehIzar2021_RNA.h5ad']

# metadata_df = metadata_df[metadata_df['Dataset'] != 'McFarlandTsherniak2020.h5ad']

metadata_df = metadata_df.sort_values(['Size (GB)'])

import sys
sys.path.append("/Genomics/pritykinlab/dillon/software/slurm_submitter")
import slurm_submitter

slurm_params = {
    'cpus-per-task': 1,
    'mem-per-cpu': '32GB',
    # Add more parameters as needed
}

# Prepare the list of arguments for each task
arguments_list = [{'row':row, 'default_slurm_params':slurm_params} for idx, row in metadata_df.iterrows()]


results = slurm_submitter.run(run_parallel.run, arguments_list, slurm_params=slurm_params, run_type='regular')

