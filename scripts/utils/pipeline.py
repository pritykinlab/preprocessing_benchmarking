import os
import scanpy as sc
import pandas as pd
import numpy as np
import scipy
from scipy.sparse import csr_matrix
from itertools import product
from sklearn.neighbors import NearestNeighbors
import sys
sys.path.append("/Genomics/pritykinlab/dillon/software/slurm_submitter")
import slurm_submitter
import itertools
import time

def run_pipeline(input_adata_file, output_dir, pipeline_params, default_slurm_params, verbose=False, parallel_type="slurm"):
    """  Executes Pipeline on scRNA-seq data based on provided steps and parameters.
    input_adata_file:
        Adata file must have a raw input
    
    """

    aggregated_filename = os.path.join(output_dir, "aggregated_results.tsv")
    if os.path.exists(aggregated_filename):
        print("Already Created aggregated data, so skipping run")
        return

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    if verbose:
        print(f"Running pipeline with {len(pipeline_params)} steps")
    intermediate_files_output_dir = os.path.join(output_dir, "intermediate_files")
    if not os.path.exists(intermediate_files_output_dir):
        os.makedirs(intermediate_files_output_dir)
    cleaned_input_adata_file = prepare_cleaned_input(input_adata_file, intermediate_files_output_dir)
    if verbose:
        print(f"Cleaned input file Completed")
        print(f"Intermediate files will be saved in {intermediate_files_output_dir}")
    
    ##################################################
    ############# First Step of Pipeline #############
    ##################################################

    all_output_adata_files_in_prev_step = [cleaned_input_adata_file]
    for i, (func, params) in enumerate(pipeline_params[:-1]):  # Exclude last step for separate processing
        print("###########################################################")
        print(f"########## Running step {i+1}/{len(pipeline_params)-1} ###############")
        print("###########################################################")

        # If file has already been written don't need to write it again assume that the result is correct
        combinations_not_complete = []
        all_combinations = itertools.product(all_output_adata_files_in_prev_step, generate_combinations(params))
        all_output_adata_files_in_prev_step = []
        for input_adata_file, comb_params in all_combinations:
            base_file_name = os.path.basename(input_adata_file).replace('.h5ad', '')
            output_file_name = f"{base_file_name}__" + "__".join([f"{key}%{value}" for key, value in comb_params.items()]) + ".h5ad"
            output_adata_file = os.path.join(intermediate_files_output_dir, output_file_name)

            if not os.path.exists(output_adata_file):
                combinations_not_complete.append([input_adata_file, comb_params])
            all_output_adata_files_in_prev_step.append(output_adata_file)

        if verbose:
            print(f"Total_combinations: {len(all_output_adata_files_in_prev_step)}")
            print(f"Combinations not complete: {combinations_not_complete}")

        # If there are combinations not complete, submit them to slurm
        arguments_l = []
        for input_adata_file, comb_params in combinations_not_complete:
            base_file_name = os.path.basename(input_adata_file).replace('.h5ad', '')
            output_file_name = f"{base_file_name}__" + "__".join([f"{key}%{value}" for key, value in comb_params.items()]) + ".h5ad"
            output_adata_file = os.path.join(intermediate_files_output_dir, output_file_name)
            arguments_l.append({"input_adata_file": input_adata_file, "output_adata_file": output_adata_file, **comb_params})

        if verbose:
            print(f"Running slurm submitter!")
        slurm_submitter.run(func, arguments_l, slurm_params=default_slurm_params, run_type=parallel_type)
        time.sleep(10) # for some reasons files take a while to write back

    ########################################################
    ############# Evaluate Section of Pipeline #############
    ########################################################
    print("###########################################################")
    print(f"########## Running Evaluation ############################")
    print("###########################################################")
    # Process the last step separately and save individual .tsv files
    all_metrics_dfs = []
    combinations_not_complete = []
    all_combinations = [(adata_file, comb_params) for adata_file in all_output_adata_files_in_prev_step for comb_params in generate_combinations(pipeline_params[-1][1])]
    print("all_combinations is", all_combinations)
    all_evaluate_files = []

    # Check which combinations need to be completed
    for adata_file, comb_params in all_combinations:
        # base_file_name = os.path.basename(adata_file).replace('.h5ad', '')
        # individual_tsv_filename = f"{base_file_name}__" + "__".join([f"{key}%{value}" for key, value in comb_params.items()]) + ".tsv"
        individual_tsv_filename = construct_individual_filename(adata_file, comb_params) + ".tsv"
        individual_tsv_path = os.path.join(intermediate_files_output_dir, individual_tsv_filename)
        if not os.path.exists(individual_tsv_path):
            combinations_not_complete.append((adata_file, comb_params))
        all_evaluate_files.append(individual_tsv_path)

    print("combinations_not_complete is", combinations_not_complete)
    print("all_evaluate_files is", all_evaluate_files)
    
    arguments_l = []
    for adata_file, comb_params in combinations_not_complete:
        # base_file_name = os.path.basename(input_adata_file).replace('.h5ad', '')
        # individual_tsv_filename = f"{base_file_name}__" + "__".join([f"{key}%{value}" for key, value in comb_params.items()]) + ".tsv"
        individual_tsv_filename = construct_individual_filename(adata_file, comb_params) + ".tsv"
        individual_tsv_path = os.path.join(intermediate_files_output_dir, individual_tsv_filename)
        arguments_l.append({"input_adata_file": adata_file, "output_file": individual_tsv_path, **comb_params})

    print("arguments_l is", arguments_l)
    
    if verbose:
        print(f"Running slurm submitter for the evaluate step!")

    func, _ = pipeline_params[-1]
    slurm_submitter.run(func, arguments_l, slurm_params=default_slurm_params, run_type=parallel_type)

    time.sleep(10) # for some reasons files take a while to write back

    # Handling the results
    for evaluate_file in all_evaluate_files:
        if os.path.exists(evaluate_file):
            metric_df = pd.read_csv(evaluate_file, sep="\t")
            params_from_file = extract_params_from_filename(evaluate_file)
            for key, value in params_from_file.items():
                metric_df[key] = value
            all_metrics_dfs.append(metric_df)
        else:
            if verbose:
                assert False, f"File {evaluate_file} does not exist!"


    # Aggregate all individual .tsv files into "aggregated_results.tsv"
    aggregated_df = pd.concat(all_metrics_dfs, ignore_index=True)
    aggregated_filename = os.path.join(output_dir, "aggregated_results.tsv")
    aggregated_df.to_csv(aggregated_filename, sep="\t", index=False)
    os.system("rm -r {intermediate_files_output_dir}")
    return aggregated_df

def construct_individual_filename(adata_file, comb_params):
    """
    Construct a unique filename for the individual .tsv files based on the input parameters.
    """
    base_name = os.path.basename(adata_file).replace('.h5ad', '')
    param_str = "__".join([f"{key}%{value}" for key, value in comb_params.items()])
    return f"{base_name}__{param_str}"


def generate_combinations(params):
    """Generate all combinations of parameters."""
    keys, values = zip(*params.items())
    return [dict(zip(keys, v)) for v in product(*values)]

def extract_params_from_filename(filename):
    """Extract parameters and their values from a filename."""
    # Initialize an empty dictionary to store extracted parameters
    params = {}

    # Remove the file extension and any directory paths
    base_filename = os.path.splitext(os.path.basename(filename))[0]

    # Split the base filename into key-value pairs on the double underscore '__'
    kv_pairs = base_filename.split('__')

    # For each key-value pair, split on the percent '%' to separate the key from the value
    for kv in kv_pairs:
        if '%' in kv:
            key, value = kv.split('%', 1)  # Only split on the first percent
            params[key] = value

    return params

def prepare_cleaned_input(input_adata_file, output_dir):
    """Prepare 'cleaned_input.h5ad' from the input AnnData object."""
    output_adata_file = os.path.join(output_dir, "cleaned_input.h5ad")
    if not os.path.exists(output_adata_file):
        adata = sc.read_h5ad(input_adata_file)
        if not isinstance(adata.X, scipy.sparse.spmatrix):
            adata.X = csr_matrix(adata.X)
        adata.layers['raw'] = adata.X
        sc.pp.filter_cells(adata, min_genes=100)
        adata.write_h5ad(output_adata_file)
        print("Finished preparing 'cleaned_input.h5ad'")
    return output_adata_file



