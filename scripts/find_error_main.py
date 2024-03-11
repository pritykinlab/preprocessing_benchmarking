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

output_base_dir = "/Genomics/pritykinlab/yujie/preprocessing_benchmarking/results"
input_datasets_dir = "/Genomics/pritykinlab/yujie/preprocessing_benchmarking/datasets/harmonized_perturb_datasets"

metadata_df = pd.read_csv("../dataset_metadata/selected_datasets_summary.tsv", sep='\t')
metadata_df = metadata_df.sort_values(['Size (GB)'])

import sys
sys.path.append("/Genomics/pritykinlab/dillon/software/slurm_submitter")
import slurm_submitter



# Prepare the list of arguments for each task
arguments_list = [{'row':row} for idx, row in metadata_df.iterrows()]

slurm_params = {
    'cpus-per-task': 2,
    'mem-per-cpu': '64GB',
    # Add more parameters as needed
}

results = slurm_submitter.run(run_parallel.run, arguments_list, slurm_params=slurm_params, run_type='slurm')

