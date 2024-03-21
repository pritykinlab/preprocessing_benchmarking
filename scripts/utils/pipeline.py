import os
import scanpy as sc
import pandas as pd
import numpy as np
import scipy
from itertools import product
import sys
sys.path.append("/Genomics/pritykinlab/dillon/software/slurm_submitter")
import slurm_submitter
import itertools
import time
from scipy.sparse import csr_matrix



import itertools
import os


def run_pipeline(input_adata_file, output_dir, pipeline_params, default_slurm_params, verbose=False, parallel_type="slurm", remove_intermediate=True):
    """  Executes Pipeline on scRNA-seq data based on provided steps and parameters.
    input_adata_file:
        Adata file must have a raw input
    """
    for pipeline_step, params in pipeline_params:
        assert callable(pipeline_step), f"Step {pipeline_step} is not callable"
        assert isinstance(params, dict), f"Params {params} is not a dictionary"
        for key, value in params.items():
            if not isinstance(value, list):
                params[key] = [value]

    original_input_adata_file = input_adata_file
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


    # Path for intermediate files
    intermediate_dir = os.path.join(output_dir, "intermediate_files")
    
    # Ensure the intermediate files directory exists
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)

    # Initialize a list to store previous combinations
    prev_combinations = [()]

    # if mapping file is already done then we don't need to get the mapping files

    mapping_files_completed = True
    for step_num in range(len(pipeline_params)):
        mapping_filename = os.path.join(intermediate_dir, f"{step_num + 1}_mapping.txt")
        if not os.path.exists(mapping_filename):
            mapping_files_completed = False
            break
    
    if not mapping_files_completed:
        print("Creating Mapping Files")
        # Process each step
        for step_num, (step, params) in enumerate(pipeline_params):
            # Generate all combinations of parameters for this step
            step_combinations = list(itertools.product(prev_combinations, *params.values()))
            
            # Flatten and combine with previous combinations
            combined_combinations = []
            for prev_combo, *current_combo in step_combinations:
                combined_combinations.append(prev_combo + tuple(current_combo))

            # Update prev_combinations for next iteration
            prev_combinations = combined_combinations

            # Determine the file extension for the current step
            file_extension = ".tsv" if step_num == len(pipeline_params) - 1 else ".h5ad"

            # Create a mapping file for this step
            mapping_filename = os.path.join(intermediate_dir, f"{step_num + 1}_mapping.txt")
            with open(mapping_filename, 'w') as file:
                # Determine all parameter names up to this step
                param_names = [name for _, step_params in pipeline_params[:step_num+1] for name in step_params]
                file.write("file_name\t" + "\t".join(param_names) + "\n")  # Header

                # Write each combination of parameters and the associated filename
                for combo_num, combo in enumerate(combined_combinations, start=1):
                    file_name = f"{step_num + 1}_{combo_num}{file_extension}"
                    param_values = "\t".join(map(str, combo))
                    file.write(f"{file_name}\t{param_values}\n")
    else:
        print("Mapping Files already created")
    

    print("#####################################################")
    print("#################### On Execution #######################")
    print("#####################################################")

    # Execution part using the mapping files
    for step_num in range(len(pipeline_params)):
        print("#############################################################")
        print(f"#################### stemp_num {step_num} #######################")
        print("##############################################################")

        step_func, current_params = pipeline_params[step_num]
        mapping_filename = os.path.join(intermediate_dir, f"{step_num + 1}_mapping.txt")

        df = pd.read_csv(mapping_filename, sep='\t')
        if step_num > 0:
            previous_df = pd.read_csv(os.path.join(intermediate_dir, f"{step_num}_mapping.txt"), sep='\t')
        else:
            previous_df = None

        arguments_l = []
        for _, params in df.iterrows():
            output_adata_file = os.path.join(intermediate_dir, params['file_name'])
            if os.path.exists(output_adata_file):
                print(f"File {output_adata_file} already exists, skipping")
                continue

            # Determine the correct input file
            if step_num == 0:
                input_adata_file = original_input_adata_file  # Use the initial input file for the first step
            else:
                # Extract previous parameters (all but those in the current step)
                previous_params_keys = df.columns.difference(current_params.keys()).tolist()
                previous_params_keys.remove('file_name')
                previous_params = params[previous_params_keys]
                matching_row = previous_df[(previous_df[previous_params_keys] == previous_params.values).all(axis=1)]
                input_adata_file = os.path.join(intermediate_dir, matching_row.iloc[0]['file_name'])

            # Extract only the parameters relevant for the current combination
            comb_params = params[current_params.keys()].to_dict()

            arguments_l.append({
                "input_adata_file": input_adata_file,
                "output_file": output_adata_file,
                **comb_params
            })


        # Run the step if there are combinations to process
        if verbose:
            print(f"Running step {step_num + 1} with {len(arguments_l)} combinations")
        if arguments_l:
            slurm_submitter.run(step_func, arguments_l, slurm_params=default_slurm_params, run_type=parallel_type)
        time.sleep(30)

    # Handling the results
    all_metrics_dfs = []
    last_mapping_filename = os.path.join(intermediate_dir, f"{len(pipeline_params)}_mapping.txt")

    # Read the last mapping file as a DataFrame
    last_mapping_df = pd.read_csv(last_mapping_filename, sep='\t')
    for _, row in last_mapping_df.iterrows():
        evaluate_file = os.path.join(intermediate_dir, row['file_name'])
        if os.path.exists(evaluate_file):
            metric_df = pd.read_csv(evaluate_file, sep="\t")
            
            # Extract all the parameters from that row, excluding 'file_name'
            params_from_file = row.drop('file_name').to_dict()
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

    if remove_intermediate:
        os.system(f"rm -r {intermediate_files_output_dir}")
    return aggregated_df





