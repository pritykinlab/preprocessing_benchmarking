from utils import pipeline
from utils import processing_steps
import os
import scipy
from itertools import product
from utils import run_parallel
import glob


output_base_dir = "/Genomics/pritykinlab/yujie/preprocessing_benchmarking/10x_CITE_seq"
input_datasets_dir = "/Genomics/pritykinlab/yujie/preprocessing_benchmarking/datasets/CITE_seq/10x_CITE_seq/h5ad"

import sys
sys.path.append("/Genomics/pritykinlab/dillon/software/slurm_submitter")
import slurm_submitter


slurm_params = {
    'cpus-per-task': 1,
    'mem-per-cpu': '32GB',
    # Add more parameters as needed
}

# Prepare the list of arguments for each task
arguments_list = [{'dataset':dataset, 'default_slurm_params':slurm_params} 
                    for dataset in glob.glob(input_datasets_dir + "/*.h5ad")]

print(arguments_list)

error_log_file = "/Genomics/pritykinlab/dillon/preprocessing_benchmarking/scripts/error.log"

results = slurm_submitter.run(run_parallel.run_10x_CITE_seq,
                              arguments_list,
                              slurm_params=slurm_params,
                              run_type='slurm',
                              max_num_jobs=5,
                              error_log_file=error_log_file)


